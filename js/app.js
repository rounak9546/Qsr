/**
 * app.js - Master Application Orchestrator & Pipeline Lifecycle Manager
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const App = (function () {
  'use strict';

  // Application state
  let currentModel = null;
  let currentScaler = null;
  let currentModelType = 'linear';
  let trainedFeatures = [];
  let modelBenchmarkHistory = [];
  let lastTestActual = [];
  let lastTestPredicted = [];
  let lastTestResiduals = [];
  let rankedCandidates = [];

  // Default features available
  const ALL_FEATURES = [
    { key: 'molecular_weight', label: 'Molecular Weight (MW)' },
    { key: 'logP', label: 'Calculated LogP' },
    { key: 'hbd', label: 'H-Bond Donors (HBD)' },
    { key: 'hba', label: 'H-Bond Acceptors (HBA)' },
    { key: 'rotatable_bonds', label: 'Rotatable Bonds' },
    { key: 'tpsa', label: 'Polar Surface Area (TPSA)' },
    { key: 'ring_count', label: 'Ring Count' },
    { key: 'aromatic_ring_count', label: 'Aromatic Rings' },
    { key: 'fraction_sp3', label: 'Fraction Csp3' }
  ];

  /**
   * Sample compound presets for rapid testing and academic demonstration
   */
  const SAMPLE_COMPOUNDS = [
    {
      name: 'Gefitinib (Iressa)',
      smiles: 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1',
      desc: '1st Gen EGFR TKI approved for non-small cell lung cancer (NSCLC)'
    },
    {
      name: 'Erlotinib (Tarceva)',
      smiles: 'COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC',
      desc: '1st Gen EGFR inhibitor with high affinity for ATP binding cleft'
    },
    {
      name: 'Osimertinib (Tagrisso)',
      smiles: 'CN(C)CC=CC(=O)Nc1cc(Nc2ncccn2)c(OC)cc1Nc3ccc(C)c(N(C)C)c3',
      desc: '3rd Gen mutant-selective EGFR inhibitor targeting T790M resistance'
    },
    {
      name: 'Lapatinib (Tykerb)',
      smiles: 'CS(=O)(=O)CCNCc1ccc(-o1)c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2',
      desc: 'Dual EGFR / HER2 receptor tyrosine kinase inhibitor'
    },
    {
      name: 'Afatinib (Gilotrif)',
      smiles: 'CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(Cl)c(F)c3)ncnc2cc1OC4CCOC4',
      desc: '2nd Gen irreversible covalent ErbB family blocker'
    },
    {
      name: 'Novel 4-Anilinoquinazoline Candidate (Hypothetical)',
      smiles: 'COc1cc2c(Nc3ccc(F)c(Cl)c3)c(C#N)cnc2cc1OCCCN4CCN(C)CC4',
      desc: 'Design variant with 3-cyanoquinoline core and basic side-chain'
    }
  ];

  /**
   * Initialize Application
   */
  async function init() {
    console.log('Initializing EGFR QSAR Bioinformatics Platform...');
    setupNavigation();
    setupModals();
    setupCSVUpload();
    setupDatasetControls();
    setupQSARTrainingControls();
    setupPredictionControls();
    setupCandidatePrioritization();
    setupExportControls();
    setupVivaAccordion();

    // Check initial URL hash
    const initialHash = window.location.hash ? window.location.hash.substring(1) : 'overview';
    UI.navigateTo(initialHash);

    // Auto-load benchmark dataset
    try {
      await loadBenchmarkDataset();
      // Pre-train default baseline model for immediate interactive demo
      trainDefaultBaselineModel();
    } catch (err) {
      console.warn('Initial benchmark dataset fetch encountered an issue:', err);
      UI.showToast('Ready for CSV upload or click "Load Benchmark Dataset".', 'info');
    }
  }

  /**
   * Set up SPA Navigation
   */
  function setupNavigation() {
    const navItems = document.querySelectorAll('.nav-item');
    navItems.forEach(item => {
      item.addEventListener('click', e => {
        e.preventDefault();
        const target = item.getAttribute('data-target');
        if (target) {
          UI.navigateTo(target);
        }
      });
    });

    // Mobile hamburger toggle
    const hamburger = document.getElementById('mobile-menu-btn');
    const sidebar = document.getElementById('app-sidebar');
    if (hamburger && sidebar) {
      hamburger.addEventListener('click', () => {
        sidebar.classList.toggle('mobile-open');
      });
    }

    // Quick workflow action buttons on overview dashboard
    const workflowButtons = document.querySelectorAll('[data-workflow-nav]');
    workflowButtons.forEach(btn => {
      btn.addEventListener('click', () => {
        const dest = btn.getAttribute('data-workflow-nav');
        if (dest) UI.navigateTo(dest);
      });
    });
  }

  /**
   * Set up Compound Inspection Modal
   */
  function setupModals() {
    const modal = document.getElementById('compound-modal');
    const closeBtn = document.getElementById('modal-close-btn');
    if (modal && closeBtn) {
      closeBtn.addEventListener('click', () => {
        UI.closeModal();
      });
      modal.addEventListener('click', e => {
        if (e.target === modal) {
          UI.closeModal();
        }
      });
    }
  }

  /**
   * Load Default Benchmark Dataset
   */
  async function loadBenchmarkDataset() {
    UI.showToast('Loading EGFR benchmark dataset from ChEMBL...', 'info', 2500);
    try {
      const result = await DataManager.loadDefaultDataset();
      onDatasetLoaded(result.dataset, result.stats);
      UI.showToast(`Loaded ${result.dataset.length} EGFR inhibitor compounds successfully.`, 'success');
    } catch (err) {
      console.error(err);
      UI.showToast('Could not load local dataset file. Please upload a CSV.', 'warning');
    }
  }

  /**
   * CSV Upload and Drag & Drop handler
   */
  function setupCSVUpload() {
    const dropZone = document.getElementById('csv-drop-zone');
    const fileInput = document.getElementById('csv-file-input');
    const loadDefaultBtn = document.getElementById('load-default-dataset-btn');

    if (loadDefaultBtn) {
      loadDefaultBtn.addEventListener('click', () => {
        loadBenchmarkDataset();
      });
    }

    if (!dropZone || !fileInput) return;

    dropZone.addEventListener('click', () => fileInput.click());

    dropZone.addEventListener('dragover', e => {
      e.preventDefault();
      dropZone.classList.add('drop-active');
    });

    dropZone.addEventListener('dragleave', () => {
      dropZone.classList.remove('drop-active');
    });

    dropZone.addEventListener('drop', e => {
      e.preventDefault();
      dropZone.classList.remove('drop-active');
      if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
        handleUploadedFile(e.dataTransfer.files[0]);
      }
    });

    fileInput.addEventListener('change', e => {
      if (e.target.files && e.target.files.length > 0) {
        handleUploadedFile(e.target.files[0]);
      }
    });
  }

  function handleUploadedFile(file) {
    if (!file.name.endsWith('.csv')) {
      UI.showToast('Please upload a valid CSV file (.csv).', 'error');
      return;
    }

    const reader = new FileReader();
    reader.onload = async function (e) {
      try {
        const text = e.target.result;
        if (!text || text.trim().length === 0) {
          UI.showToast('Uploaded CSV file is empty.', 'error');
          return;
        }

        const rawData = await DataManager.parseCSV(text);
        if (!rawData || rawData.length === 0) {
          UI.showToast('Failed to parse CSV records. Please check file formatting.', 'error');
          return;
        }

        // Validate headers
        const headers = Object.keys(rawData[0]);
        const validation = Validator.validateCSVHeaders(headers);
        if (!validation.valid) {
          UI.showToast(`CSV missing required columns: ${validation.missing.join(', ')}`, 'error', 6000);
          return;
        }

        const result = DataManager.preprocessDataset(rawData);
        if (result.dataset.length === 0) {
          UI.showToast('No valid compound records could be extracted from this CSV.', 'error');
          return;
        }

        onDatasetLoaded(result.dataset, result.stats);
        UI.showToast(`Successfully processed ${result.dataset.length} compounds from uploaded CSV.`, 'success');
      } catch (err) {
        console.error(err);
        UI.showToast(`CSV Processing Error: ${err.message}`, 'error');
      }
    };
    reader.onerror = function () {
      UI.showToast('Error reading the uploaded file from disk.', 'error');
    };
    reader.readAsText(file);
  }

  /**
   * Action taken when a dataset is successfully loaded or refreshed
   */
  function onDatasetLoaded(dataset, stats) {
    // Update Overview stats cards
    updateOverviewStats(dataset);

    // Update Preprocessing Summary Panel
    updatePreprocessingDisplay(stats);

    // Render Table in Dataset Explorer
    UI.getDatasetTableState().page = 1;
    UI.renderDatasetTable(dataset);

    // Populate Molecular Analysis Charts
    updateMolecularAnalysis(dataset);

    // Reset Model state when dataset changes
    currentModel = null;
    currentScaler = null;
    resetModelResultsView();
  }

  /**
   * Update Overview Dashboard KPI cards
   */
  function updateOverviewStats(dataset) {
    const compCountEl = document.getElementById('kpi-compounds-count');
    const descCountEl = document.getElementById('kpi-descriptors-count');
    const activeCountEl = document.getElementById('kpi-active-count');
    const targetEl = document.getElementById('kpi-target-label');

    if (compCountEl) compCountEl.textContent = dataset.length;
    if (descCountEl) descCountEl.textContent = '9 Descriptors';
    if (targetEl) targetEl.textContent = 'EGFR (ErbB1 / HER1)';

    if (activeCountEl) {
      const activeCount = dataset.filter(d => d.pic50 >= 7.0).length;
      activeCountEl.textContent = `${activeCount} Potent Leads`;
    }
  }

  /**
   * Update Preprocessing Display Card
   */
  function updatePreprocessingDisplay(stats) {
    const rawCountEl = document.getElementById('stat-raw-rows');
    const dupCountEl = document.getElementById('stat-duplicates-removed');
    const invalidCountEl = document.getElementById('stat-invalid-removed');
    const finalCountEl = document.getElementById('stat-final-rows');
    const logListEl = document.getElementById('preprocessing-log-list');

    if (rawCountEl) rawCountEl.textContent = stats.originalRows;
    if (dupCountEl) dupCountEl.textContent = stats.duplicatesRemoved;
    if (invalidCountEl) invalidCountEl.textContent = stats.invalidRows;
    if (finalCountEl) finalCountEl.textContent = stats.finalRows;

    if (logListEl) {
      logListEl.innerHTML = stats.log.map(msg => `<li>${Validator.sanitize(msg)}</li>`).join('');
    }
  }

  /**
   * Setup controls for Dataset Explorer table (Search, Sort, Filter, Pagination)
   */
  function setupDatasetControls() {
    const searchInput = document.getElementById('dataset-search-input');
    const filterSelect = document.getElementById('dataset-activity-filter');
    const pageSizeSelect = document.getElementById('dataset-page-size');
    const prevBtn = document.getElementById('dataset-prev-page');
    const nextBtn = document.getElementById('dataset-next-page');

    if (searchInput) {
      searchInput.addEventListener('input', e => {
        UI.getDatasetTableState().searchQuery = e.target.value;
        UI.getDatasetTableState().page = 1;
        UI.renderDatasetTable(DataManager.getCurrentDataset());
      });
    }

    if (filterSelect) {
      filterSelect.addEventListener('change', e => {
        UI.getDatasetTableState().activityFilter = e.target.value;
        UI.getDatasetTableState().page = 1;
        UI.renderDatasetTable(DataManager.getCurrentDataset());
      });
    }

    if (pageSizeSelect) {
      pageSizeSelect.addEventListener('change', e => {
        UI.getDatasetTableState().pageSize = parseInt(e.target.value, 10) || 10;
        UI.getDatasetTableState().page = 1;
        UI.renderDatasetTable(DataManager.getCurrentDataset());
      });
    }

    if (prevBtn) {
      prevBtn.addEventListener('click', () => {
        if (UI.getDatasetTableState().page > 1) {
          UI.getDatasetTableState().page--;
          UI.renderDatasetTable(DataManager.getCurrentDataset());
        }
      });
    }

    if (nextBtn) {
      nextBtn.addEventListener('click', () => {
        UI.getDatasetTableState().page++;
        UI.renderDatasetTable(DataManager.getCurrentDataset());
      });
    }

    // Column sorting headers
    const sortHeaders = document.querySelectorAll('#dataset-table-head th[data-sort]');
    sortHeaders.forEach(th => {
      th.addEventListener('click', () => {
        const col = th.getAttribute('data-sort');
        const state = UI.getDatasetTableState();
        if (state.sortColumn === col) {
          state.sortDirection = state.sortDirection === 'asc' ? 'desc' : 'asc';
        } else {
          state.sortColumn = col;
          state.sortDirection = 'desc';
        }
        UI.renderDatasetTable(DataManager.getCurrentDataset());
      });
    });
  }

  /**
   * Update Molecular Analysis Visualizations
   */
  function updateMolecularAnalysis(dataset) {
    if (!dataset || dataset.length === 0) return;

    // 1. Activity pIC50 distribution
    const pic50s = dataset.map(d => d.pic50);
    Charts.renderDistribution('chart-pic50-dist', pic50s, 'Biological Activity (pIC₅₀)', '-log₁₀ M', 8);

    // 2. Molecular Weight distribution
    const mws = dataset.map(d => d.molecular_weight);
    Charts.renderDistribution('chart-mw-dist', mws, 'Molecular Weight', 'Da', 8);

    // 3. LogP distribution
    const logPs = dataset.map(d => d.logP);
    Charts.renderDistribution('chart-logp-dist', logPs, 'Calculated LogP', '', 8);

    // 4. TPSA distribution
    const tpsas = dataset.map(d => d.tpsa);
    Charts.renderDistribution('chart-tpsa-dist', tpsas, 'Polar Surface Area (TPSA)', 'Å²', 8);

    // 5. PCA Projection
    const featureKeys = ALL_FEATURES.map(f => f.key);
    const pcaPoints = DataManager.computePCA(dataset, featureKeys);
    Charts.renderPCAPlot('chart-pca', pcaPoints);

    // 6. Correlation Heatmap Table
    const corr = DataManager.computeCorrelationMatrix(dataset, featureKeys);
    renderCorrelationHeatmapTable(corr);
  }

  /**
   * Render Correlation Matrix as an interactive grid
   */
  function renderCorrelationHeatmapTable(corr) {
    const container = document.getElementById('correlation-matrix-container');
    if (!container) return;

    const shortLabels = ['MW', 'LogP', 'HBD', 'HBA', 'RotB', 'TPSA', 'Rings', 'Arom', 'Fsp3'];

    let html = `<table class="corr-table"><thead><tr><th></th>`;
    shortLabels.forEach(lbl => {
      html += `<th>${lbl}</th>`;
    });
    html += `</tr></thead><tbody>`;

    for (let i = 0; i < corr.matrix.length; i++) {
      html += `<tr><th>${shortLabels[i]}</th>`;
      for (let j = 0; j < corr.matrix[i].length; j++) {
        const val = corr.matrix[i][j];
        // Color mapping: red for negative, cyan/teal for positive
        let bg = 'transparent';
        if (val > 0) {
          bg = `rgba(6, 182, 212, ${Math.min(0.85, Math.abs(val))})`;
        } else if (val < 0) {
          bg = `rgba(239, 68, 68, ${Math.min(0.85, Math.abs(val))})`;
        }
        html += `<td style="background-color: ${bg};" title="${corr.features[i]} vs ${corr.features[j]}: r = ${val}">${val.toFixed(2)}</td>`;
      }
      html += `</tr>`;
    }
    html += `</tbody></table>`;
    container.innerHTML = html;
  }

  /**
   * Set up QSAR Model Training Form
   */
  function setupQSARTrainingControls() {
    const trainBtn = document.getElementById('train-model-btn');
    const featureContainer = document.getElementById('feature-checkbox-container');
    const selectAllBtn = document.getElementById('select-all-features-btn');
    const deselectAllBtn = document.getElementById('deselect-all-features-btn');

    // Populate feature checkboxes
    if (featureContainer) {
      featureContainer.innerHTML = ALL_FEATURES.map(
        f => `
        <label class="checkbox-label">
          <input type="checkbox" name="qsar-feature" value="${f.key}" checked />
          <span>${f.label}</span>
        </label>
      `
      ).join('');
    }

    if (selectAllBtn) {
      selectAllBtn.addEventListener('click', () => {
        const boxes = document.querySelectorAll('input[name="qsar-feature"]');
        boxes.forEach(b => (b.checked = true));
      });
    }

    if (deselectAllBtn) {
      deselectAllBtn.addEventListener('click', () => {
        const boxes = document.querySelectorAll('input[name="qsar-feature"]');
        boxes.forEach(b => (b.checked = false));
      });
    }

    if (trainBtn) {
      trainBtn.addEventListener('click', () => {
        executeModelTraining();
      });
    }
  }

  /**
   * Train Default Baseline Model on App Initialization
   */
  function trainDefaultBaselineModel() {
    const dataset = DataManager.getCurrentDataset();
    if (!dataset || dataset.length < 10) return;

    console.log('Training default baseline QSAR model (Random Forest & Linear Regression)...');
    executeModelTraining(false); // silent mode
  }

  /**
   * Execute QSAR Model Training Workflow
   */
  function executeModelTraining(notify = true) {
    const dataset = DataManager.getCurrentDataset();
    if (!dataset || dataset.length === 0) {
      UI.showToast('Cannot train model: No compound dataset loaded.', 'error');
      return;
    }

    // Selected features
    const checkedBoxes = document.querySelectorAll('input[name="qsar-feature"]:checked');
    const selectedFeatures = Array.from(checkedBoxes).map(b => b.value);

    if (selectedFeatures.length === 0) {
      UI.showToast('Please select at least 1 molecular descriptor as an input feature.', 'warning');
      return;
    }

    // Hyperparameters from UI
    const modelSelect = document.getElementById('model-algorithm-select');
    const splitSelect = document.getElementById('train-test-split-select');
    const seedInput = document.getElementById('random-seed-input');
    const cvSelect = document.getElementById('cross-validation-select');

    const algorithm = modelSelect ? modelSelect.value : 'randomForest';
    const splitRatio = splitSelect ? parseFloat(splitSelect.value) : 0.2;
    const seed = seedInput ? parseInt(seedInput.value, 10) || 42 : 42;
    const cvFolds = cvSelect ? parseInt(cvSelect.value, 10) || 0 : 5;

    currentModelType = algorithm;
    trainedFeatures = selectedFeatures;

    // Build feature matrix X and target y (pIC50)
    const X = dataset.map(row => selectedFeatures.map(feat => Number(row[feat]) || 0));
    const y = dataset.map(row => Number(row.pic50));

    // Train/Test Split
    const split = ML.trainTestSplit(X, y, splitRatio, seed);

    // Standardize features using training data only to avoid data leakage
    const scaler = new ML.StandardScaler();
    const XTrainScaled = scaler.fitTransform(split.XTrain);
    const XTestScaled = scaler.transform(split.XTest);
    currentScaler = scaler;

    let modelInstance = null;
    let modelDisplayName = 'Multiple Linear Regression';

    if (algorithm === 'linear') {
      modelInstance = new ML.LinearRegressionModel(1e-4);
      modelDisplayName = 'Multiple Linear Regression (OLS/Ridge)';
    } else if (algorithm === 'ridge') {
      modelInstance = new ML.LinearRegressionModel(1e-2);
      modelDisplayName = 'Ridge Regression (L2 = 0.01)';
    } else if (algorithm === 'randomForest') {
      modelInstance = new ML.RandomForestRegressorModel(20, 6, 3, 0.75, seed);
      modelDisplayName = 'Random Forest Regressor (20 Trees)';
    }

    // Fit model
    const featureLabels = selectedFeatures.map(f => {
      const match = ALL_FEATURES.find(af => af.key === f);
      return match ? match.label : f;
    });

    modelInstance.fit(XTrainScaled, split.yTrain, featureLabels);
    currentModel = modelInstance;

    // Predictions on Training and Testing sets
    const yPredTrain = modelInstance.predict(XTrainScaled);
    const yPredTest = modelInstance.predict(XTestScaled);

    const trainMetrics = ML.calculateRegressionMetrics(split.yTrain, yPredTrain);
    const testMetrics = ML.calculateRegressionMetrics(split.yTest, yPredTest);

    lastTestActual = split.yTest;
    lastTestPredicted = yPredTest;
    lastTestResiduals = testMetrics.residuals;

    // Perform K-Fold Cross Validation if selected
    let cvResults = null;
    if (cvFolds > 0) {
      try {
        cvResults = ML.kFoldCrossValidation(
          algorithm === 'randomForest' ? 'randomForest' : 'linear',
          X,
          y,
          cvFolds,
          { seed, featureNames: featureLabels }
        );
      } catch (cvErr) {
        console.warn('Cross-validation warning:', cvErr);
      }
    }

    // Add to benchmark history
    const historyItem = {
      id: Date.now(),
      name: modelDisplayName,
      algorithm: algorithm,
      trainR2: trainMetrics.r2,
      r2: testMetrics.r2,
      rmse: testMetrics.rmse,
      mae: testMetrics.mae,
      featuresCount: selectedFeatures.length,
      trainSize: split.XTrain.length,
      testSize: split.XTest.length,
      cvMean: cvResults ? cvResults.meanScore : null,
      cvStd: cvResults ? cvResults.stdScore : null,
      cvFolds: cvFolds
    };

    modelBenchmarkHistory.unshift(historyItem);

    // Update UI Results Dashboard
    updateModelResultsDashboard(historyItem, modelInstance, testMetrics, trainMetrics, cvResults);

    // Automatically run prioritization on current dataset
    prioritizeCandidates(false);

    if (notify) {
      UI.showToast(`Model "${modelDisplayName}" trained successfully! Test R² = ${testMetrics.r2.toFixed(3)}`, 'success');
      UI.navigateTo('model-results');
    }
  }

  /**
   * Update Model Results Dashboard with KPIs, Charts, and Comparison Table
   */
  function updateModelResultsDashboard(modelSummary, model, testMetrics, trainMetrics, cvResults) {
    // KPI Cards
    const r2El = document.getElementById('metric-r2-val');
    const rmseEl = document.getElementById('metric-rmse-val');
    const maeEl = document.getElementById('metric-mae-val');
    const mseEl = document.getElementById('metric-mse-val');
    const trainR2El = document.getElementById('metric-train-r2-val');
    const splitInfoEl = document.getElementById('model-split-info');
    const cvInfoEl = document.getElementById('model-cv-info');

    if (r2El) r2El.textContent = testMetrics.r2.toFixed(3);
    if (rmseEl) rmseEl.textContent = testMetrics.rmse.toFixed(3);
    if (maeEl) maeEl.textContent = testMetrics.mae.toFixed(3);
    if (mseEl) mseEl.textContent = testMetrics.mse.toFixed(3);
    if (trainR2El) trainR2El.textContent = trainMetrics.r2.toFixed(3);

    if (splitInfoEl) {
      splitInfoEl.textContent = `Training: ${modelSummary.trainSize} compounds | Test: ${modelSummary.testSize} compounds`;
    }

    if (cvInfoEl) {
      if (cvResults) {
        cvInfoEl.innerHTML = `<strong>${cvResults.k}-Fold Cross-Validation:</strong> Mean R² = ${cvResults.meanScore.toFixed(3)} (± ${cvResults.stdScore.toFixed(3)})`;
      } else {
        cvInfoEl.textContent = 'Cross-validation was disabled for this run.';
      }
    }

    // Hide empty state banner
    const emptyState = document.getElementById('model-results-empty-state');
    const resultsContainer = document.getElementById('model-results-active-container');
    if (emptyState) emptyState.style.display = 'none';
    if (resultsContainer) resultsContainer.style.display = 'block';

    // 1. Render Actual vs Predicted Chart
    Charts.renderActualVsPredicted('chart-actual-vs-pred', lastTestActual, lastTestPredicted, modelSummary.name);

    // 2. Render Residual Plot
    Charts.renderResidualPlot('chart-residuals', lastTestPredicted, lastTestResiduals, modelSummary.name);

    // 3. Render Feature Importance Bar Chart
    const importance = model.getFeatureImportance();
    Charts.renderFeatureImportance('chart-feature-importance', importance);

    // 4. Update Model Comparison Table and Benchmark Chart
    updateModelComparison();
  }

  /**
   * Reset Results View when dataset is reset or empty
   */
  function resetModelResultsView() {
    const emptyState = document.getElementById('model-results-empty-state');
    const resultsContainer = document.getElementById('model-results-active-container');
    if (emptyState) emptyState.style.display = 'block';
    if (resultsContainer) resultsContainer.style.display = 'none';
  }

  /**
   * Update Model Benchmark Comparison Table and Chart
   */
  function updateModelComparison() {
    const tableBody = document.getElementById('model-comparison-table-body');
    if (!tableBody) return;

    if (modelBenchmarkHistory.length === 0) {
      tableBody.innerHTML = `<tr><td colspan="7" class="text-center py-4 text-muted">No models evaluated yet.</td></tr>`;
      return;
    }

    // Find best model by highest test R2
    let bestR2 = -Infinity;
    let bestId = null;
    modelBenchmarkHistory.forEach(m => {
      if (m.r2 > bestR2) {
        bestR2 = m.r2;
        bestId = m.id;
      }
    });

    tableBody.innerHTML = modelBenchmarkHistory
      .map(m => {
        const isBest = m.id === bestId;
        const cvDisplay = m.cvMean !== null ? `${m.cvMean.toFixed(3)} (±${m.cvStd.toFixed(3)})` : '—';

        return `
        <tr class="${isBest ? 'highlight-best-row' : ''}">
          <td><strong>${Validator.sanitize(m.name)}</strong> ${isBest ? '<span class="badge badge-active ml-2">Highest Test R²</span>' : ''}</td>
          <td><strong>${m.r2.toFixed(3)}</strong></td>
          <td>${m.rmse.toFixed(3)}</td>
          <td>${m.mae.toFixed(3)}</td>
          <td>${m.trainR2.toFixed(3)}</td>
          <td>${cvDisplay}</td>
          <td>${m.trainSize}/${m.testSize}</td>
        </tr>
      `;
      })
      .join('');

    // Render benchmark chart
    Charts.renderModelComparisonChart('chart-model-benchmark', modelBenchmarkHistory);
  }

  /**
   * Set up Interactive Prediction Section Controls
   */
  function setupPredictionControls() {
    const smilesInput = document.getElementById('predict-smiles-input');
    const nameInput = document.getElementById('predict-name-input');
    const predictBtn = document.getElementById('run-prediction-btn');
    const sampleSelect = document.getElementById('predict-sample-select');

    // Populate preset dropdown
    if (sampleSelect) {
      sampleSelect.innerHTML =
        `<option value="">-- Choose an authentic EGFR compound or lead --</option>` +
        SAMPLE_COMPOUNDS.map(
          (c, idx) => `<option value="${idx}">${Validator.sanitize(c.name)}</option>`
        ).join('');

      sampleSelect.addEventListener('change', e => {
        const idx = e.target.value;
        if (idx !== '' && SAMPLE_COMPOUNDS[idx]) {
          const comp = SAMPLE_COMPOUNDS[idx];
          if (nameInput) nameInput.value = comp.name;
          if (smilesInput) smilesInput.value = comp.smiles;
          // Live calculate and render structure
          handleSingleCompoundPrediction();
        }
      });
    }

    if (predictBtn) {
      predictBtn.addEventListener('click', () => {
        handleSingleCompoundPrediction();
      });
    }

    if (smilesInput) {
      smilesInput.addEventListener('keypress', e => {
        if (e.key === 'Enter') {
          handleSingleCompoundPrediction();
        }
      });
    }
  }

  /**
   * Handle Single Compound Activity Prediction
   */
  function handleSingleCompoundPrediction() {
    const smilesInput = document.getElementById('predict-smiles-input');
    const nameInput = document.getElementById('predict-name-input');
    const outputContainer = document.getElementById('prediction-output-card');
    const canvas = document.getElementById('predict-structure-canvas');

    if (!smilesInput) return;
    const smiles = smilesInput.value.trim();
    const name = nameInput && nameInput.value.trim() ? nameInput.value.trim() : 'Query Compound';

    if (!smiles) {
      UI.showToast('Please enter a SMILES molecular string.', 'warning');
      return;
    }

    // Validate SMILES
    const validation = Validator.validateSMILES(smiles);
    if (!validation.valid) {
      UI.showToast(`Invalid SMILES structure: ${validation.error}`, 'error', 5000);
      return;
    }

    // Calculate Descriptors
    const desc = Descriptors.calculateDescriptorsFromSMILES(smiles);

    // Render 2D Chemical Structure
    if (canvas) {
      Descriptors.renderMolecule(canvas, smiles, { bgColor: '#0b1120' });
    }

    // Check if model exists
    if (!currentModel || !currentScaler) {
      UI.showToast('Please train a QSAR model in the "QSAR Modeling" section first.', 'warning');
      if (outputContainer) {
        outputContainer.innerHTML = `
          <div class="alert alert-warning">
            <strong>Model Not Trained:</strong> Please navigate to the QSAR Modeling tab and click "Train QSAR Model" to enable predictive activity scoring.
          </div>
        `;
        outputContainer.style.display = 'block';
      }
      return;
    }

    // Prepare feature vector aligned with trainedFeatures
    const featureVector = trainedFeatures.map(featKey => desc[featKey] || 0);

    // Scale feature vector
    const scaledVector = currentScaler.transformRow(featureVector);

    // Model Prediction
    const predictedPIC50 = currentModel.predictRow(scaledVector);
    const estIC50 = QSAR.pic50ToNanomolar(predictedPIC50);
    const classification = QSAR.classifyActivity(predictedPIC50);

    // Evaluate Applicability Domain
    const rawTrainData = DataManager.getCurrentDataset().map(row =>
      trainedFeatures.map(k => Number(row[k]) || 0)
    );
    const scaledTrainData = currentScaler.transform(rawTrainData);
    const adEval = QSAR.evaluateApplicabilityDomain(scaledTrainData, scaledVector);

    // Render Prediction Output Card
    if (outputContainer) {
      outputContainer.innerHTML = `
        <div class="prediction-result-header">
          <div>
            <h3 class="text-xl font-bold text-slate-100">${Validator.sanitize(name)}</h3>
            <span class="font-mono text-xs text-muted break-all">${Validator.sanitize(smiles)}</span>
          </div>
          <span class="badge ${classification.badgeClass}">${classification.label}</span>
        </div>

        <div class="prediction-metrics-banner my-4">
          <div class="metric-block">
            <span class="metric-label">Model-Predicted pIC₅₀</span>
            <span class="metric-val text-cyan">${predictedPIC50.toFixed(2)}</span>
            <small class="text-muted">-log₁₀(Molar IC₅₀)</small>
          </div>
          <div class="metric-block">
            <span class="metric-label">Estimated Potency (IC₅₀)</span>
            <span class="metric-val text-teal">${estIC50 < 1000 ? estIC50.toFixed(1) + ' nM' : (estIC50 / 1000).toFixed(2) + ' µM'}</span>
            <small class="text-muted">Calculated from pIC₅₀</small>
          </div>
          <div class="metric-block">
            <span class="metric-label">Applicability Domain</span>
            <span class="metric-val ${adEval.insideDomain ? 'text-emerald' : 'text-amber'}">${adEval.insideDomain ? 'Inside' : 'Borderline / Outside'}</span>
            <small class="text-muted">Distance: ${adEval.distance.toFixed(2)} (Limit: ${adEval.threshold.toFixed(2)})</small>
          </div>
        </div>

        <div class="applicability-domain-notice ${adEval.insideDomain ? 'ad-pass' : 'ad-warn'} p-3 rounded mb-4">
          <p class="text-xs">${Validator.sanitize(adEval.message)}</p>
        </div>

        <h4 class="text-sm font-semibold text-cyan uppercase tracking-wider mb-2">Calculated Molecular Descriptors</h4>
        <div class="descriptors-summary-table mb-4">
          <div class="desc-item"><span>Molecular Weight:</span> <strong>${desc.molecular_weight.toFixed(2)} Da</strong></div>
          <div class="desc-item"><span>LogP:</span> <strong>${desc.logP.toFixed(2)}</strong></div>
          <div class="desc-item"><span>H-Bond Donors:</span> <strong>${desc.hbd}</strong></div>
          <div class="desc-item"><span>H-Bond Acceptors:</span> <strong>${desc.hba}</strong></div>
          <div class="desc-item"><span>Rotatable Bonds:</span> <strong>${desc.rotatable_bonds}</strong></div>
          <div class="desc-item"><span>TPSA:</span> <strong>${desc.tpsa.toFixed(2)} Å²</strong></div>
          <div class="desc-item"><span>Ring Count:</span> <strong>${desc.ring_count}</strong></div>
          <div class="desc-item"><span>Aromatic Rings:</span> <strong>${desc.aromatic_ring_count}</strong></div>
          <div class="desc-item"><span>Fraction Csp3:</span> <strong>${desc.fraction_sp3.toFixed(2)}</strong></div>
        </div>

        <div class="lipinski-box ${desc.lipinski_violations <= 1 ? 'lipinski-pass' : 'lipinski-fail'} mb-3">
          <strong>Lipinski Rule of 5:</strong> ${desc.lipinski_violations === 0 ? '✓ 0 Violations (Good oral druglikeness potential)' : `⚠️ ${desc.lipinski_violations} Violation(s) [${desc.lipinski_details.join(', ')}]`}
        </div>

        <div class="scientific-disclaimer-box">
          <span class="disclaimer-badge">Academic & Scientific Notice</span>
          <p class="text-xs text-muted mt-1">
            This activity value is a <strong>computationally prioritized model estimation</strong> trained on EGFR kinase bioassays. Computational predictions do not prove clinical anticancer efficacy, cellular selectivity, or pharmacokinetic safety. Formal biochemical kinase inhibition assays (e.g. TR-FRET, ADP-Glo) and in vitro cell viability validations are strictly required.
          </p>
        </div>
      `;
      outputContainer.style.display = 'block';
    }

    UI.showToast(`Predicted pIC₅₀ = ${predictedPIC50.toFixed(2)} for ${name}`, 'success');
  }

  /**
   * Set up Candidate Prioritization / Virtual Screening Controls
   */
  function setupCandidatePrioritization() {
    const rankBtn = document.getElementById('run-candidate-ranking-btn');
    const searchInput = document.getElementById('candidate-search-input');
    const pageSizeSelect = document.getElementById('candidate-page-size');
    const prevBtn = document.getElementById('candidate-prev-page');
    const nextBtn = document.getElementById('candidate-next-page');

    if (rankBtn) {
      rankBtn.addEventListener('click', () => {
        prioritizeCandidates(true);
      });
    }

    if (searchInput) {
      searchInput.addEventListener('input', e => {
        UI.getCandidateTableState().searchQuery = e.target.value;
        UI.getCandidateTableState().page = 1;
        UI.renderCandidateTable(rankedCandidates);
      });
    }

    if (pageSizeSelect) {
      pageSizeSelect.addEventListener('change', e => {
        UI.getCandidateTableState().pageSize = parseInt(e.target.value, 10) || 10;
        UI.getCandidateTableState().page = 1;
        UI.renderCandidateTable(rankedCandidates);
      });
    }

    if (prevBtn) {
      prevBtn.addEventListener('click', () => {
        if (UI.getCandidateTableState().page > 1) {
          UI.getCandidateTableState().page--;
          UI.renderCandidateTable(rankedCandidates);
        }
      });
    }

    if (nextBtn) {
      nextBtn.addEventListener('click', () => {
        UI.getCandidateTableState().page++;
        UI.renderCandidateTable(rankedCandidates);
      });
    }

    // Sorting headers for candidate table
    const sortHeaders = document.querySelectorAll('#candidate-table-head th[data-sort]');
    sortHeaders.forEach(th => {
      th.addEventListener('click', () => {
        const col = th.getAttribute('data-sort');
        const state = UI.getCandidateTableState();
        if (state.sortColumn === col) {
          state.sortDirection = state.sortDirection === 'asc' ? 'desc' : 'asc';
        } else {
          state.sortColumn = col;
          state.sortDirection = 'desc';
        }
        UI.renderCandidateTable(rankedCandidates);
      });
    });
  }

  /**
   * Run Candidate Prioritization across the entire dataset
   */
  function prioritizeCandidates(notify = true) {
    const dataset = DataManager.getCurrentDataset();
    if (!dataset || dataset.length === 0) {
      if (notify) UI.showToast('No dataset available to prioritize.', 'error');
      return;
    }

    if (!currentModel || !currentScaler) {
      if (notify) UI.showToast('Please train a QSAR model first before ranking candidates.', 'warning');
      return;
    }

    const ranked = dataset.map(row => {
      const fVec = trainedFeatures.map(k => Number(row[k]) || 0);
      const scaled = currentScaler.transformRow(fVec);
      const pred = currentModel.predictRow(scaled);

      return {
        compound_id: row.compound_id,
        compound_name: row.compound_name,
        smiles: row.smiles,
        actual_pic50: row.pic50,
        predicted_pic50: parseFloat(pred.toFixed(2)),
        activity_class: row.activity_class
      };
    });

    // Sort descending by predicted pIC50
    ranked.sort((a, b) => b.predicted_pic50 - a.predicted_pic50);
    rankedCandidates = ranked;

    UI.getCandidateTableState().page = 1;
    UI.renderCandidateTable(rankedCandidates);

    if (notify) {
      UI.showToast(`Prioritized and ranked ${ranked.length} compounds by predicted activity.`, 'success');
      UI.navigateTo('candidate-ranking');
    }
  }

  /**
   * Set up Export buttons
   */
  function setupExportControls() {
    const exportDatasetBtn = document.getElementById('export-dataset-csv-btn');
    const exportCandidatesBtn = document.getElementById('export-candidates-csv-btn');
    const exportMetricsBtn = document.getElementById('export-metrics-csv-btn');

    if (exportDatasetBtn) {
      exportDatasetBtn.addEventListener('click', () => {
        const dataset = DataManager.getCurrentDataset();
        if (!dataset || dataset.length === 0) {
          UI.showToast('No dataset to export.', 'warning');
          return;
        }
        DataManager.exportToCSV('egfr_processed_dataset.csv', dataset);
        UI.showToast('Processed dataset exported to CSV.', 'success');
      });
    }

    if (exportCandidatesBtn) {
      exportCandidatesBtn.addEventListener('click', () => {
        if (!rankedCandidates || rankedCandidates.length === 0) {
          UI.showToast('No prioritized candidates to export. Run candidate ranking first.', 'warning');
          return;
        }
        DataManager.exportToCSV('egfr_prioritized_candidates.csv', rankedCandidates);
        UI.showToast('Prioritized candidates exported to CSV.', 'success');
      });
    }

    if (exportMetricsBtn) {
      exportMetricsBtn.addEventListener('click', () => {
        if (!modelBenchmarkHistory || modelBenchmarkHistory.length === 0) {
          UI.showToast('No model metrics to export. Train at least one model.', 'warning');
          return;
        }
        DataManager.exportToCSV('egfr_qsar_model_metrics.csv', modelBenchmarkHistory);
        UI.showToast('Model metrics exported to CSV.', 'success');
      });
    }
  }

  /**
   * Set up Viva Voce FAQ Accordion & Search Filter
   */
  function setupVivaAccordion() {
    const searchInput = document.getElementById('viva-search-input');
    const questions = document.querySelectorAll('.viva-card');

    questions.forEach(card => {
      const header = card.querySelector('.viva-header');
      if (header) {
        header.addEventListener('click', () => {
          card.classList.toggle('open');
        });
      }
    });

    if (searchInput) {
      searchInput.addEventListener('input', e => {
        const query = e.target.value.toLowerCase().trim();
        questions.forEach(card => {
          const text = card.textContent.toLowerCase();
          if (!query || text.includes(query)) {
            card.style.display = 'block';
          } else {
            card.style.display = 'none';
          }
        });
      });
    }
  }

  return {
    init,
    loadBenchmarkDataset,
    executeModelTraining,
    prioritizeCandidates
  };
})();

// Bootstrap application when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
  App.init();
});
