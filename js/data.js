/**
 * data.js - Dataset Management, Preprocessing Pipeline, CSV Parser & PCA
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const DataManager = (function () {
  'use strict';

  // Active dataset state
  let currentDataset = [];
  let preprocessingStats = {
    originalRows: 0,
    duplicatesRemoved: 0,
    invalidRows: 0,
    missingValues: 0,
    finalRows: 0,
    log: []
  };

  /**
   * Robust Native CSV Parser (zero-dependency fallback if PapaParse is blocked)
   * Handles quoted commas, multiline cells, and Windows/Unix line endings
   */
  function parseCSVNative(text) {
    const lines = [];
    let row = [''];
    let inQuotes = false;

    for (let i = 0; i < text.length; i++) {
      const char = text[i];
      const nextChar = text[i + 1];

      if (char === '"') {
        if (inQuotes && nextChar === '"') {
          row[row.length - 1] += '"';
          i++; // Skip escaped quote
        } else {
          inQuotes = !inQuotes;
        }
      } else if (char === ',' && !inQuotes) {
        row.push('');
      } else if ((char === '\r' || char === '\n') && !inQuotes) {
        if (char === '\r' && nextChar === '\n') i++;
        if (row.length > 1 || row[0] !== '') {
          lines.push(row);
        }
        row = [''];
      } else {
        row[row.length - 1] += char;
      }
    }
    if (row.length > 1 || row[0] !== '') {
      lines.push(row);
    }

    if (lines.length < 2) return [];

    const rawHeaders = lines[0].map(h => h.trim());
    const data = [];

    for (let i = 1; i < lines.length; i++) {
      const line = lines[i];
      if (line.length === 0 || (line.length === 1 && line[0].trim() === '')) continue;
      const obj = {};
      for (let j = 0; j < rawHeaders.length; j++) {
        const header = rawHeaders[j];
        obj[header] = line[j] !== undefined ? line[j].trim() : '';
      }
      data.push(obj);
    }

    return data;
  }

  /**
   * Load CSV string using PapaParse if present, otherwise native fallback
   */
  function parseCSV(csvString) {
    return new Promise((resolve, reject) => {
      try {
        if (typeof window.Papa !== 'undefined') {
          window.Papa.parse(csvString, {
            header: true,
            skipEmptyLines: true,
            dynamicTyping: false,
            complete: function (results) {
              resolve(results.data);
            },
            error: function (err) {
              console.warn('PapaParse error, using native parser fallback:', err);
              resolve(parseCSVNative(csvString));
            }
          });
        } else {
          resolve(parseCSVNative(csvString));
        }
      } catch (e) {
        try {
          resolve(parseCSVNative(csvString));
        } catch (inner) {
          reject(inner);
        }
      }
    });
  }

  /**
   * Preprocess and clean raw dataset
   * Performs:
   * 1. Header canonicalization
   * 2. Duplicate detection (by SMILES and compound_id)
   * 3. Missing/invalid values detection
   * 4. Activity unit conversion to Molar and pIC50 calculation
   * 5. Descriptor calculation if missing
   */
  function preprocessDataset(rawData) {
    const stats = {
      originalRows: rawData.length,
      duplicatesRemoved: 0,
      invalidRows: 0,
      missingValues: 0,
      finalRows: 0,
      log: []
    };

    stats.log.push(`Ingested raw dataset with ${rawData.length} records.`);

    const cleaned = [];
    const seenSmiles = new Set();
    const seenIds = new Set();

    rawData.forEach((row, index) => {
      // Find case-insensitive keys
      const keys = Object.keys(row);
      const findKey = patterns => keys.find(k => patterns.some(p => k.toLowerCase().replace(/[\s_-]/g, '').includes(p)));

      const smilesKey = findKey(['smiles']);
      const idKey = findKey(['compoundid', 'chemblid', 'id']);
      const nameKey = findKey(['compoundname', 'name', 'pref_name']);
      const actValKey = findKey(['activityvalue', 'standardvalue', 'ic50', 'value']);
      const actUnitKey = findKey(['activityunit', 'standardunits', 'unit']);
      const pic50Key = findKey(['pic50']);

      const smiles = smilesKey ? String(row[smilesKey]).trim() : '';
      const compoundId = idKey ? String(row[idKey]).trim() : `CMPD-${index + 1}`;
      const compoundName = nameKey ? String(row[nameKey]).trim() : compoundId;

      // 1. SMILES check
      if (!smiles) {
        stats.invalidRows++;
        return;
      }

      // 2. Duplicate check
      if (seenSmiles.has(smiles) || (compoundId && seenIds.has(compoundId))) {
        stats.duplicatesRemoved++;
        return;
      }
      seenSmiles.add(smiles);
      if (compoundId) seenIds.add(compoundId);

      // 3. Activity value & pIC50 derivation
      let pic50 = null;
      let actVal = actValKey ? parseFloat(row[actValKey]) : NaN;
      let actUnit = actUnitKey ? String(row[actUnitKey]).trim() : 'nM';

      if (pic50Key && !isNaN(parseFloat(row[pic50Key]))) {
        pic50 = parseFloat(row[pic50Key]);
      } else if (!isNaN(actVal) && actVal > 0) {
        const conv = QSAR.calculatePIC50(actVal, actUnit);
        if (conv) {
          pic50 = conv.pic50;
        }
      }

      if (pic50 === null || isNaN(pic50) || pic50 <= 0 || pic50 > 14) {
        stats.invalidRows++;
        return;
      }

      // 4. Extract or compute Descriptors
      const mwKey = findKey(['molecularweight', 'mw', 'molwt']);
      const logPKey = findKey(['logp', 'alogp']);
      const hbdKey = findKey(['hbd', 'hdonors']);
      const hbaKey = findKey(['hba', 'hacceptors']);
      const rotKey = findKey(['rotatablebonds', 'rotbonds']);
      const tpsaKey = findKey(['tpsa', 'polararea']);
      const ringKey = findKey(['ringcount', 'rings']);
      const aromKey = findKey(['aromaticringcount', 'aromaticrings']);
      const sp3Key = findKey(['fractionssp3', 'fsp3', 'fractions_sp3']);

      let mw = mwKey && !isNaN(parseFloat(row[mwKey])) ? parseFloat(row[mwKey]) : null;
      let logP = logPKey && !isNaN(parseFloat(row[logPKey])) ? parseFloat(row[logPKey]) : null;
      let hbd = hbdKey && !isNaN(parseInt(row[hbdKey], 10)) ? parseInt(row[hbdKey], 10) : null;
      let hba = hbaKey && !isNaN(parseInt(row[hbaKey], 10)) ? parseInt(row[hbaKey], 10) : null;
      let rot = rotKey && !isNaN(parseInt(row[rotKey], 10)) ? parseInt(row[rotKey], 10) : null;
      let tpsa = tpsaKey && !isNaN(parseFloat(row[tpsaKey])) ? parseFloat(row[tpsaKey]) : null;
      let rings = ringKey && !isNaN(parseInt(row[ringKey], 10)) ? parseInt(row[ringKey], 10) : null;
      let aromRings = aromKey && !isNaN(parseInt(row[aromKey], 10)) ? parseInt(row[aromKey], 10) : null;
      let fsp3 = sp3Key && !isNaN(parseFloat(row[sp3Key])) ? parseFloat(row[sp3Key]) : null;

      // If any descriptor is missing, compute dynamically from SMILES
      if (mw === null || logP === null || hbd === null || hba === null || tpsa === null) {
        const computed = Descriptors.calculateDescriptorsFromSMILES(smiles);
        if (mw === null) mw = computed.molecular_weight;
        if (logP === null) logP = computed.logP;
        if (hbd === null) hbd = computed.hbd;
        if (hba === null) hba = computed.hba;
        if (rot === null) rot = computed.rotatable_bonds;
        if (tpsa === null) tpsa = computed.tpsa;
        if (rings === null) rings = computed.ring_count;
        if (aromRings === null) aromRings = computed.aromatic_ring_count;
        if (fsp3 === null) fsp3 = computed.fraction_sp3;
      }

      // If IC50 value wasn't provided, estimate from pIC50
      if (isNaN(actVal) || actVal <= 0) {
        actVal = QSAR.pic50ToNanomolar(pic50);
        actUnit = 'nM';
      }

      const activityClassification = QSAR.classifyActivity(pic50);

      cleaned.push({
        compound_id: compoundId,
        compound_name: compoundName,
        smiles: smiles,
        target: 'EGFR',
        activity_type: 'IC50',
        activity_value: parseFloat(actVal.toFixed(2)),
        activity_unit: actUnit,
        pic50: parseFloat(pic50.toFixed(2)),
        activity_class: activityClassification.label,
        binary_class: activityClassification.binaryClass,
        molecular_weight: parseFloat(mw.toFixed(2)),
        logP: parseFloat(logP.toFixed(2)),
        hbd: hbd,
        hba: hba,
        rotatable_bonds: rot,
        tpsa: parseFloat(tpsa.toFixed(2)),
        ring_count: rings,
        aromatic_ring_count: aromRings,
        fraction_sp3: parseFloat(fsp3.toFixed(2))
      });
    });

    stats.finalRows = cleaned.length;
    stats.log.push(`Removed ${stats.duplicatesRemoved} duplicate compounds.`);
    stats.log.push(`Filtered ${stats.invalidRows} invalid or unconvertible records.`);
    stats.log.push(`Preprocessing complete. Usable analytical dataset: ${stats.finalRows} compounds.`);

    preprocessingStats = stats;
    currentDataset = cleaned;

    return {
      dataset: cleaned,
      stats: stats
    };
  }

  /**
   * Load Default Bundled EGFR Dataset from server/static path
   */
  async function loadDefaultDataset() {
    try {
      const response = await fetch('data/egfr_activity.csv');
      if (!response.ok) {
        throw new Error(`Failed to load egfr_activity.csv (Status: ${response.status})`);
      }
      const text = await response.text();
      const rawData = await parseCSV(text);
      return preprocessDataset(rawData);
    } catch (err) {
      console.error('Error fetching default dataset:', err);
      throw err;
    }
  }

  /**
   * Perform Principal Component Analysis (PCA) on continuous molecular descriptors
   * Computes top 2 Principal Components (PC1 and PC2) for 2D chemical space visualization
   * @param {Array<Object>} data
   * @param {Array<string>} featureKeys
   */
  function computePCA(data, featureKeys) {
    if (!data || data.length < 3 || !featureKeys || featureKeys.length < 2) {
      return [];
    }

    const N = data.length;
    const P = featureKeys.length;

    // Build and standardize feature matrix X
    const X = data.map(item => featureKeys.map(k => Number(item[k]) || 0));
    const scaler = new ML.StandardScaler();
    const XScaled = scaler.fitTransform(X);

    // Compute Covariance Matrix (P x P) = (X^T * X) / (N - 1)
    const cov = Array.from({ length: P }, () => new Array(P).fill(0));
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < P; j++) {
        for (let k = 0; k < P; k++) {
          cov[j][k] += (XScaled[i][j] * XScaled[i][k]) / (N - 1);
        }
      }
    }

    // Power Iteration to extract largest eigenvector (PC1)
    function powerIteration(matrix, iterations = 50) {
      let vec = new Array(P).fill(0).map(() => Math.random() - 0.5);
      for (let iter = 0; iter < iterations; iter++) {
        const nextVec = new Array(P).fill(0);
        for (let r = 0; r < P; r++) {
          for (let c = 0; c < P; c++) {
            nextVec[r] += matrix[r][c] * vec[c];
          }
        }
        // Normalize
        let norm = Math.sqrt(nextVec.reduce((sum, v) => sum + v * v, 0)) || 1.0;
        vec = nextVec.map(v => v / norm);
      }
      return vec;
    }

    const pc1Vector = powerIteration(cov);

    // Deflate covariance matrix to find PC2 orthogonal component
    const deflated = Array.from({ length: P }, () => new Array(P).fill(0));
    for (let r = 0; r < P; r++) {
      for (let c = 0; c < P; c++) {
        deflated[r][c] = cov[r][c] - pc1Vector[r] * pc1Vector[c];
      }
    }

    const pc2Vector = powerIteration(deflated);

    // Project data points onto PC1 and PC2
    return data.map((item, i) => {
      let pc1 = 0;
      let pc2 = 0;
      for (let j = 0; j < P; j++) {
        pc1 += XScaled[i][j] * pc1Vector[j];
        pc2 += XScaled[i][j] * pc2Vector[j];
      }
      return {
        id: item.compound_id,
        name: item.compound_name,
        pc1: parseFloat(pc1.toFixed(3)),
        pc2: parseFloat(pc2.toFixed(3)),
        pic50: item.pic50
      };
    });
  }

  /**
   * Compute Pearson Correlation Coefficient Matrix for feature pairs
   */
  function computeCorrelationMatrix(data, featureKeys) {
    const P = featureKeys.length;
    const N = data.length;
    const matrix = Array.from({ length: P }, () => new Array(P).fill(1));

    // Calculate means and stds
    const means = [];
    const stds = [];

    for (let j = 0; j < P; j++) {
      const vals = data.map(d => Number(d[featureKeys[j]]) || 0);
      const mean = vals.reduce((a, b) => a + b, 0) / N;
      const variance = vals.reduce((sum, v) => sum + (v - mean) * (v - mean), 0) / N;
      means.push(mean);
      stds.push(Math.sqrt(variance) || 1.0);
    }

    for (let i = 0; i < P; i++) {
      for (let j = i; j < P; j++) {
        if (i === j) {
          matrix[i][j] = 1.0;
          continue;
        }
        let cov = 0;
        for (let k = 0; k < N; k++) {
          const valI = Number(data[k][featureKeys[i]]) || 0;
          const valJ = Number(data[k][featureKeys[j]]) || 0;
          cov += (valI - means[i]) * (valJ - means[j]);
        }
        const r = cov / (N * stds[i] * stds[j]);
        const cleanR = parseFloat(r.toFixed(3));
        matrix[i][j] = cleanR;
        matrix[j][i] = cleanR;
      }
    }

    return {
      features: featureKeys,
      matrix: matrix
    };
  }

  /**
   * Export Array of Objects as CSV file download
   */
  function exportToCSV(filename, data) {
    if (!data || data.length === 0) return false;

    const headers = Object.keys(data[0]);
    let csvContent = headers.join(',') + '\n';

    data.forEach(row => {
      const line = headers
        .map(h => {
          let val = row[h] === null || row[h] === undefined ? '' : String(row[h]);
          if (val.includes(',') || val.includes('"') || val.includes('\n')) {
            val = `"${val.replace(/"/g, '""')}"`;
          }
          return val;
        })
        .join(',');
      csvContent += line + '\n';
    });

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.setAttribute('href', url);
    link.setAttribute('download', filename);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
    return true;
  }

  return {
    parseCSV,
    preprocessDataset,
    loadDefaultDataset,
    computePCA,
    computeCorrelationMatrix,
    exportToCSV,
    getCurrentDataset: () => currentDataset,
    getPreprocessingStats: () => preprocessingStats
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.DataManager = DataManager;
}
