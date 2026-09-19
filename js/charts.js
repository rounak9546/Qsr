/**
 * charts.js - Scientific Data Visualizations using Chart.js with Resilient Canvas Fallback
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const Charts = (function () {
  'use strict';

  // Global store for active Chart.js instances to avoid canvas reuse collision
  const chartInstances = {};

  // Standard scientific dark palette
  const PALETTE = {
    bg: '#0f172a',
    cardBg: '#1e293b',
    border: 'rgba(255, 255, 255, 0.08)',
    text: '#e2e8f0',
    mutedText: '#94a3b8',
    grid: 'rgba(255, 255, 255, 0.06)',
    cyan: '#06b6d4',
    teal: '#14b8a6',
    blue: '#38bdf8',
    purple: '#a855f7',
    amber: '#f59e0b',
    red: '#ef4444',
    emerald: '#10b981'
  };

  /**
   * Safely destroy previous chart instance for a canvas ID
   */
  function destroyPrevious(id) {
    if (chartInstances[id]) {
      try {
        chartInstances[id].destroy();
      } catch (e) {
        console.warn(`Could not destroy chart instance ${id}:`, e);
      }
      delete chartInstances[id];
    }
  }

  /**
   * Check if Chart.js is loaded in the browser
   */
  function isChartJsAvailable() {
    return typeof window.Chart !== 'undefined';
  }

  /**
   * Common Chart.js options for uniform scientific styling
   */
  function getBaseOptions(titleText, xLabel, yLabel) {
    return {
      responsive: true,
      maintainAspectRatio: false,
      animation: { duration: 400 },
      plugins: {
        legend: {
          labels: {
            color: PALETTE.text,
            font: { family: 'Inter, sans-serif', size: 12 }
          }
        },
        title: {
          display: !!titleText,
          text: titleText,
          color: PALETTE.text,
          font: { family: 'Inter, sans-serif', size: 14, weight: '600' },
          padding: { bottom: 12 }
        },
        tooltip: {
          backgroundColor: '#0f172a',
          titleColor: PALETTE.cyan,
          bodyColor: PALETTE.text,
          borderColor: 'rgba(6, 182, 212, 0.3)',
          borderWidth: 1,
          padding: 10
        }
      },
      scales: {
        x: {
          title: {
            display: !!xLabel,
            text: xLabel,
            color: PALETTE.mutedText,
            font: { family: 'Inter, sans-serif', size: 12 }
          },
          ticks: { color: PALETTE.mutedText },
          grid: { color: PALETTE.grid }
        },
        y: {
          title: {
            display: !!yLabel,
            text: yLabel,
            color: PALETTE.mutedText,
            font: { family: 'Inter, sans-serif', size: 12 }
          },
          ticks: { color: PALETTE.mutedText },
          grid: { color: PALETTE.grid }
        }
      }
    };
  }

  /**
   * Fallback Canvas Drawer if Chart.js is unavailable
   */
  function renderFallback(canvas, message) {
    if (!canvas || !canvas.getContext) return;
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = PALETTE.cardBg;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = PALETTE.mutedText;
    ctx.font = '13px Inter, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(message || 'Chart unavailable (Library loading or fallback mode)', canvas.width / 2, canvas.height / 2);
  }

  /**
   * 1. Actual vs Predicted Scatter Plot
   */
  function renderActualVsPredicted(canvasId, actual, predicted, modelName = 'Model') {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    destroyPrevious(canvasId);

    if (!isChartJsAvailable()) {
      renderFallback(canvas, 'Chart.js is loading or unavailable.');
      return;
    }

    const scatterData = actual.map((act, i) => ({ x: act, y: predicted[i] }));
    const minVal = Math.floor(Math.min(...actual, ...predicted) * 0.95);
    const maxVal = Math.ceil(Math.max(...actual, ...predicted) * 1.05);

    // Diagonal reference line y = x
    const diagData = [
      { x: minVal, y: minVal },
      { x: maxVal, y: maxVal }
    ];

    const options = getBaseOptions(
      `${modelName}: Actual vs Predicted pIC₅₀`,
      'Experimental pIC₅₀ (-log₁₀ M)',
      'Model-Predicted pIC₅₀ (-log₁₀ M)'
    );

    chartInstances[canvasId] = new window.Chart(canvas, {
      type: 'scatter',
      data: {
        datasets: [
          {
            label: 'Test Compounds',
            data: scatterData,
            backgroundColor: 'rgba(6, 182, 212, 0.75)',
            borderColor: PALETTE.cyan,
            borderWidth: 1.5,
            pointRadius: 5,
            pointHoverRadius: 7
          },
          {
            label: 'Identity Line (y = x)',
            type: 'line',
            data: diagData,
            borderColor: 'rgba(148, 163, 184, 0.6)',
            borderDash: [5, 5],
            borderWidth: 1.5,
            pointRadius: 0,
            fill: false
          }
        ]
      },
      options: options
    });
  }

  /**
   * 2. Residual Plot (Residuals vs Predicted)
   */
  function renderResidualPlot(canvasId, predicted, residuals, modelName = 'Model') {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    destroyPrevious(canvasId);

    if (!isChartJsAvailable()) {
      renderFallback(canvas, 'Chart.js is loading or unavailable.');
      return;
    }

    const scatterData = predicted.map((pred, i) => ({ x: pred, y: residuals[i] }));
    const minX = Math.floor(Math.min(...predicted) * 0.95);
    const maxX = Math.ceil(Math.max(...predicted) * 1.05);

    const zeroLine = [
      { x: minX, y: 0 },
      { x: maxX, y: 0 }
    ];

    const options = getBaseOptions(
      `${modelName}: Residual Plot`,
      'Predicted pIC₅₀',
      'Residual (Actual - Predicted)'
    );

    chartInstances[canvasId] = new window.Chart(canvas, {
      type: 'scatter',
      data: {
        datasets: [
          {
            label: 'Residuals',
            data: scatterData,
            backgroundColor: 'rgba(245, 158, 11, 0.75)',
            borderColor: PALETTE.amber,
            borderWidth: 1.5,
            pointRadius: 5,
            pointHoverRadius: 7
          },
          {
            label: 'Zero Residual Line',
            type: 'line',
            data: zeroLine,
            borderColor: 'rgba(239, 68, 68, 0.8)',
            borderDash: [4, 4],
            borderWidth: 1.5,
            pointRadius: 0,
            fill: false
          }
        ]
      },
      options: options
    });
  }

  /**
   * 3. Feature Importance Bar Chart
   */
  function renderFeatureImportance(canvasId, importanceList) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    destroyPrevious(canvasId);

    if (!isChartJsAvailable()) {
      renderFallback(canvas, 'Chart.js is loading or unavailable.');
      return;
    }

    if (!importanceList || importanceList.length === 0) {
      renderFallback(canvas, 'No feature importance available.');
      return;
    }

    const labels = importanceList.map(item => item.feature);
    const dataVals = importanceList.map(item => item.importance);

    const options = getBaseOptions('QSAR Feature Importance', 'Relative Importance Score', 'Descriptor');
    options.indexAxis = 'y'; // Horizontal bar chart

    chartInstances[canvasId] = new window.Chart(canvas, {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [
          {
            label: 'Relative Weight',
            data: dataVals,
            backgroundColor: 'rgba(20, 184, 166, 0.75)',
            borderColor: PALETTE.teal,
            borderWidth: 1.2,
            borderRadius: 4
          }
        ]
      },
      options: options
    });
  }

  /**
   * 4. Distribution Histogram for Descriptors
   */
  function renderDistribution(canvasId, values, label, unit = '', binCount = 8) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    destroyPrevious(canvasId);

    if (!isChartJsAvailable() || !values || values.length === 0) {
      renderFallback(canvas, 'Distribution data unavailable.');
      return;
    }

    const min = Math.min(...values);
    const max = Math.max(...values);
    const step = (max - min) / binCount || 1;

    const bins = new Array(binCount).fill(0);
    const binLabels = [];

    for (let b = 0; b < binCount; b++) {
      const bMin = min + b * step;
      const bMax = bMin + step;
      binLabels.push(`${bMin.toFixed(1)}-${bMax.toFixed(1)}`);
    }

    values.forEach(v => {
      let bIdx = Math.floor((v - min) / step);
      if (bIdx >= binCount) bIdx = binCount - 1;
      if (bIdx < 0) bIdx = 0;
      bins[bIdx]++;
    });

    const options = getBaseOptions(`${label} Distribution`, `${label} (${unit})`, 'Frequency Count');

    chartInstances[canvasId] = new window.Chart(canvas, {
      type: 'bar',
      data: {
        labels: binLabels,
        datasets: [
          {
            label: `Compounds Count`,
            data: bins,
            backgroundColor: 'rgba(56, 189, 248, 0.7)',
            borderColor: PALETTE.blue,
            borderWidth: 1.2,
            borderRadius: 4
          }
        ]
      },
      options: options
    });
  }

  /**
   * 5. Principal Component Analysis (PCA) 2D Scatter Plot
   */
  function renderPCAPlot(canvasId, pcaPoints) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    destroyPrevious(canvasId);

    if (!isChartJsAvailable() || !pcaPoints || pcaPoints.length === 0) {
      renderFallback(canvas, 'PCA data unavailable.');
      return;
    }

    // Segregate points by activity class
    const activeData = [];
    const intermediateData = [];
    const inactiveData = [];

    pcaPoints.forEach(pt => {
      const p = { x: pt.pc1, y: pt.pc2, name: pt.name };
      if (pt.pic50 >= 7.0) activeData.push(p);
      else if (pt.pic50 >= 6.0) intermediateData.push(p);
      else inactiveData.push(p);
    });

    const options = getBaseOptions('Chemical Space PCA Projection', 'Principal Component 1', 'Principal Component 2');
    options.plugins.tooltip.callbacks = {
      label: function (context) {
        const item = context.raw;
        return `${item.name || 'Compound'}: (PC1: ${item.x.toFixed(2)}, PC2: ${item.y.toFixed(2)})`;
      }
    };

    chartInstances[canvasId] = new window.Chart(canvas, {
      type: 'scatter',
      data: {
        datasets: [
          {
            label: 'Active (pIC₅₀ ≥ 7.0)',
            data: activeData,
            backgroundColor: 'rgba(16, 185, 129, 0.8)',
            borderColor: PALETTE.emerald,
            pointRadius: 5
          },
          {
            label: 'Intermediate (6.0 ≤ pIC₅₀ < 7.0)',
            data: intermediateData,
            backgroundColor: 'rgba(245, 158, 11, 0.8)',
            borderColor: PALETTE.amber,
            pointRadius: 5
          },
          {
            label: 'Inactive (pIC₅₀ < 6.0)',
            data: inactiveData,
            backgroundColor: 'rgba(239, 68, 68, 0.8)',
            borderColor: PALETTE.red,
            pointRadius: 5
          }
        ]
      },
      options: options
    });
  }

  /**
   * 6. Model Comparison Bar Chart
   */
  function renderModelComparisonChart(canvasId, modelResults) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    destroyPrevious(canvasId);

    if (!isChartJsAvailable() || !modelResults || modelResults.length === 0) {
      renderFallback(canvas, 'Train models to generate comparative performance metrics.');
      return;
    }

    const labels = modelResults.map(m => m.name);
    const r2Vals = modelResults.map(m => m.r2);
    const rmseVals = modelResults.map(m => m.rmse);
    const maeVals = modelResults.map(m => m.mae);

    const options = getBaseOptions('Model Benchmark Comparison', 'Algorithm', 'Score / Error');

    chartInstances[canvasId] = new window.Chart(canvas, {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [
          {
            label: 'Test R²',
            data: r2Vals,
            backgroundColor: 'rgba(6, 182, 212, 0.75)',
            borderColor: PALETTE.cyan,
            borderWidth: 1.2,
            borderRadius: 4
          },
          {
            label: 'RMSE',
            data: rmseVals,
            backgroundColor: 'rgba(245, 158, 11, 0.75)',
            borderColor: PALETTE.amber,
            borderWidth: 1.2,
            borderRadius: 4
          },
          {
            label: 'MAE',
            data: maeVals,
            backgroundColor: 'rgba(168, 85, 247, 0.75)',
            borderColor: PALETTE.purple,
            borderWidth: 1.2,
            borderRadius: 4
          }
        ]
      },
      options: options
    });
  }

  return {
    destroyPrevious,
    renderActualVsPredicted,
    renderResidualPlot,
    renderFeatureImportance,
    renderDistribution,
    renderPCAPlot,
    renderModelComparisonChart,
    renderFallback,
    isChartJsAvailable
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.Charts = Charts;
}
