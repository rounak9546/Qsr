/**
 * ml.js - Machine Learning Core (Linear Regression, Ridge, Random Forest, Cross-Validation, Metrics)
 * Pure JavaScript implementation - Fast, deterministic, and 100% self-contained
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const ML = (function () {
  'use strict';

  /**
   * Mulberry32 seeded Pseudo-Random Number Generator (PRNG)
   * Ensures deterministic, reproducible train/test splits and random forest bootstrap
   */
  function createPRNG(seed = 42) {
    let s = Math.floor(seed) >>> 0;
    return function () {
      s |= 0;
      s = (s + 0x6d2b79f5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /**
   * Shuffle array in-place using seeded PRNG
   */
  function shuffleArray(arr, prng) {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(prng() * (i + 1));
      const temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
    return arr;
  }

  /**
   * Standard Scaler (Z-Score Normalization)
   * Fitted strictly on training set to prevent data leakage
   */
  class StandardScaler {
    constructor() {
      this.means = [];
      this.stds = [];
      this.isFitted = false;
    }

    fit(X) {
      const N = X.length;
      if (N === 0) return;
      const P = X[0].length;

      this.means = new Array(P).fill(0);
      this.stds = new Array(P).fill(0);

      // Compute feature means
      for (let i = 0; i < N; i++) {
        for (let j = 0; j < P; j++) {
          this.means[j] += X[i][j];
        }
      }
      for (let j = 0; j < P; j++) {
        this.means[j] /= N;
      }

      // Compute standard deviations
      for (let i = 0; i < N; i++) {
        for (let j = 0; j < P; j++) {
          const diff = X[i][j] - this.means[j];
          this.stds[j] += diff * diff;
        }
      }
      for (let j = 0; j < P; j++) {
        const variance = this.stds[j] / N;
        this.stds[j] = Math.sqrt(variance) || 1.0; // avoid division by zero
      }

      this.isFitted = true;
    }

    transform(X) {
      if (!this.isFitted) throw new Error('StandardScaler must be fitted before transform.');
      return X.map(row =>
        row.map((val, j) => (val - this.means[j]) / this.stds[j])
      );
    }

    fitTransform(X) {
      this.fit(X);
      return this.transform(X);
    }

    transformRow(row) {
      if (!this.isFitted) throw new Error('StandardScaler must be fitted.');
      return row.map((val, j) => (val - this.means[j]) / this.stds[j]);
    }
  }

  /**
   * Train/Test Split
   * @param {Array} X - features matrix
   * @param {Array} y - target array
   * @param {number} testSize - fraction between 0.1 and 0.5 (e.g. 0.2 for 80/20)
   * @param {number} seed - random seed for reproducibility
   */
  function trainTestSplit(X, y, testSize = 0.2, seed = 42) {
    if (!X || !y || X.length !== y.length) {
      throw new Error('X and y must be non-empty and of matching length.');
    }

    const n = X.length;
    const indices = Array.from({ length: n }, (_, i) => i);
    const prng = createPRNG(seed);
    shuffleArray(indices, prng);

    const testCount = Math.max(1, Math.round(n * testSize));
    const trainIndices = indices.slice(testCount);
    const testIndices = indices.slice(0, testCount);

    const XTrain = trainIndices.map(i => X[i]);
    const yTrain = trainIndices.map(i => y[i]);
    const XTest = testIndices.map(i => X[i]);
    const yTest = testIndices.map(i => y[i]);

    return {
      XTrain,
      yTrain,
      XTest,
      yTest,
      trainIndices,
      testIndices
    };
  }

  /**
   * Linear Algebra: Solve (A^T * A + lambda * I) * w = A^T * b
   * Uses Gaussian elimination with partial pivoting
   */
  function solveLinearSystem(A, b) {
    const n = A.length;
    // Augmented matrix
    const M = A.map((row, i) => [...row, b[i]]);

    for (let i = 0; i < n; i++) {
      // Find pivot
      let maxEl = Math.abs(M[i][i]);
      let maxRow = i;
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(M[k][i]) > maxEl) {
          maxEl = Math.abs(M[k][i]);
          maxRow = k;
        }
      }

      // Swap rows
      if (maxRow !== i) {
        const tmp = M[i];
        M[i] = M[maxRow];
        M[maxRow] = tmp;
      }

      if (Math.abs(M[i][i]) < 1e-12) {
        M[i][i] = 1e-12; // Numerical stability regularizer
      }

      // Eliminate below
      for (let k = i + 1; k < n; k++) {
        const c = -M[k][i] / M[i][i];
        for (let j = i; j <= n; j++) {
          if (i === j) {
            M[k][j] = 0;
          } else {
            M[k][j] += c * M[i][j];
          }
        }
      }
    }

    // Back substitution
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = M[i][n] / M[i][i];
      for (let k = i - 1; k >= 0; k--) {
        M[k][n] -= M[k][i] * x[i];
      }
    }

    return x;
  }

  /**
   * Linear / Ridge Regression Model
   */
  class LinearRegressionModel {
    constructor(alpha = 1e-4) {
      this.alpha = alpha; // L2 Regularization coefficient
      this.weights = [];
      this.bias = 0;
      this.featureNames = [];
    }

    fit(X, y, featureNames = []) {
      const N = X.length;
      const P = X[0].length;
      this.featureNames = featureNames;

      // Add intercept term to feature matrix: X_aug has shape N x (P + 1)
      // Normal equations: (X_aug^T * X_aug + alpha * I) * theta = X_aug^T * y
      const dim = P + 1;
      const XtX = Array.from({ length: dim }, () => new Array(dim).fill(0));
      const Xty = new Array(dim).fill(0);

      for (let i = 0; i < N; i++) {
        const xRow = [...X[i], 1.0]; // bias feature is 1.0
        const yi = y[i];

        for (let r = 0; r < dim; r++) {
          Xty[r] += xRow[r] * yi;
          for (let c = 0; c < dim; c++) {
            XtX[r][c] += xRow[r] * xRow[c];
          }
        }
      }

      // Add Ridge regularization (do not penalize intercept at index P)
      for (let j = 0; j < P; j++) {
        XtX[j][j] += this.alpha * N;
      }

      const theta = solveLinearSystem(XtX, Xty);
      this.weights = theta.slice(0, P);
      this.bias = theta[P];
    }

    predict(X) {
      return X.map(row => {
        let yHat = this.bias;
        for (let j = 0; j < row.length; j++) {
          yHat += this.weights[j] * row[j];
        }
        return yHat;
      });
    }

    predictRow(row) {
      let yHat = this.bias;
      for (let j = 0; j < row.length; j++) {
        yHat += this.weights[j] * row[j];
      }
      return yHat;
    }

    getFeatureImportance() {
      // For standardized features, absolute magnitude of weights corresponds to feature importance
      const total = this.weights.reduce((sum, w) => sum + Math.abs(w), 0) || 1.0;
      return this.weights.map((w, i) => ({
        feature: this.featureNames[i] || `Feature ${i + 1}`,
        weight: parseFloat(w.toFixed(4)),
        importance: parseFloat((Math.abs(w) / total).toFixed(4))
      })).sort((a, b) => b.importance - a.importance);
    }
  }

  /**
   * Decision Tree Node for Random Forest Regressor
   */
  class DecisionTreeNode {
    constructor() {
      this.isLeaf = false;
      this.value = 0;
      this.featureIndex = null;
      this.threshold = null;
      this.left = null;
      this.right = null;
    }
  }

  /**
   * Random Forest Regressor
   * Pure JavaScript implementation of Breiman's Random Forest algorithm
   */
  class RandomForestRegressorModel {
    constructor(numTrees = 20, maxDepth = 6, minSamplesSplit = 3, maxFeaturesRatio = 0.7, seed = 42) {
      this.numTrees = numTrees;
      this.maxDepth = maxDepth;
      this.minSamplesSplit = minSamplesSplit;
      this.maxFeaturesRatio = maxFeaturesRatio;
      this.seed = seed;
      this.trees = [];
      this.featureNames = [];
      this.featureImportanceScores = [];
    }

    fit(X, y, featureNames = []) {
      const N = X.length;
      const P = X[0].length;
      this.featureNames = featureNames;
      this.featureImportanceScores = new Array(P).fill(0);
      this.trees = [];

      const prng = createPRNG(this.seed);
      const mTry = Math.max(1, Math.round(P * this.maxFeaturesRatio));

      for (let t = 0; t < this.numTrees; t++) {
        // Bootstrap sample (sampling with replacement)
        const sampleIndices = [];
        for (let i = 0; i < N; i++) {
          sampleIndices.push(Math.floor(prng() * N));
        }

        const tree = this.buildTree(X, y, sampleIndices, 0, mTry, prng);
        this.trees.push(tree);
      }

      // Normalize feature importance
      const totalImportance = this.featureImportanceScores.reduce((a, b) => a + b, 0) || 1.0;
      this.normalizedImportance = this.featureImportanceScores.map((score, i) => ({
        feature: this.featureNames[i] || `Feature ${i + 1}`,
        weight: parseFloat(score.toFixed(4)),
        importance: parseFloat((score / totalImportance).toFixed(4))
      })).sort((a, b) => b.importance - a.importance);
    }

    buildTree(X, y, sampleIndices, depth, mTry, prng) {
      const node = new DecisionTreeNode();
      const nSamples = sampleIndices.length;

      // Base cases: max depth reached, insufficient samples, or zero variance
      let sum = 0;
      for (let i = 0; i < nSamples; i++) {
        sum += y[sampleIndices[i]];
      }
      const meanVal = sum / nSamples;

      if (depth >= this.maxDepth || nSamples < this.minSamplesSplit) {
        node.isLeaf = true;
        node.value = meanVal;
        return node;
      }

      // Calculate current node variance (MSE * N)
      let parentVarianceSum = 0;
      for (let i = 0; i < nSamples; i++) {
        const diff = y[sampleIndices[i]] - meanVal;
        parentVarianceSum += diff * diff;
      }

      if (parentVarianceSum < 1e-7) {
        node.isLeaf = true;
        node.value = meanVal;
        return node;
      }

      // Select random subset of features
      const P = X[0].length;
      const allFeatures = Array.from({ length: P }, (_, i) => i);
      shuffleArray(allFeatures, prng);
      const candidateFeatures = allFeatures.slice(0, mTry);

      let bestGain = -Infinity;
      let bestFeature = null;
      let bestThreshold = null;
      let bestLeftIndices = null;
      let bestRightIndices = null;

      for (let f of candidateFeatures) {
        // Extract values for candidate feature
        const featVals = sampleIndices.map(idx => X[idx][f]);
        featVals.sort((a, b) => a - b);

        // Test potential threshold splits
        const step = Math.max(1, Math.floor(nSamples / 10));
        for (let k = 1; k < nSamples; k += step) {
          if (featVals[k] === featVals[k - 1]) continue;
          const threshold = (featVals[k] + featVals[k - 1]) / 2;

          const leftIndices = [];
          const rightIndices = [];
          let leftSum = 0;
          let rightSum = 0;

          for (let i = 0; i < nSamples; i++) {
            const idx = sampleIndices[i];
            const val = X[idx][f];
            const targetVal = y[idx];
            if (val <= threshold) {
              leftIndices.push(idx);
              leftSum += targetVal;
            } else {
              rightIndices.push(idx);
              rightSum += targetVal;
            }
          }

          if (leftIndices.length === 0 || rightIndices.length === 0) continue;

          const leftMean = leftSum / leftIndices.length;
          const rightMean = rightSum / rightIndices.length;

          let leftVarianceSum = 0;
          for (let idx of leftIndices) {
            const d = y[idx] - leftMean;
            leftVarianceSum += d * d;
          }

          let rightVarianceSum = 0;
          for (let idx of rightIndices) {
            const d = y[idx] - rightMean;
            rightVarianceSum += d * d;
          }

          const gain = parentVarianceSum - (leftVarianceSum + rightVarianceSum);

          if (gain > bestGain) {
            bestGain = gain;
            bestFeature = f;
            bestThreshold = threshold;
            bestLeftIndices = leftIndices;
            bestRightIndices = rightIndices;
          }
        }
      }

      if (bestGain <= 0 || !bestLeftIndices || !bestRightIndices) {
        node.isLeaf = true;
        node.value = meanVal;
        return node;
      }

      // Record feature importance gain
      this.featureImportanceScores[bestFeature] += bestGain;

      node.featureIndex = bestFeature;
      node.threshold = bestThreshold;
      node.left = this.buildTree(X, y, bestLeftIndices, depth + 1, mTry, prng);
      node.right = this.buildTree(X, y, bestRightIndices, depth + 1, mTry, prng);

      return node;
    }

    predictSample(tree, row) {
      let curr = tree;
      while (!curr.isLeaf) {
        if (row[curr.featureIndex] <= curr.threshold) {
          curr = curr.left;
        } else {
          curr = curr.right;
        }
      }
      return curr.value;
    }

    predict(X) {
      return X.map(row => this.predictRow(row));
    }

    predictRow(row) {
      let sum = 0;
      for (let t = 0; t < this.trees.length; t++) {
        sum += this.predictSample(this.trees[t], row);
      }
      return sum / this.trees.length;
    }

    getFeatureImportance() {
      return this.normalizedImportance || [];
    }
  }

  /**
   * Logistic Regression Model (Binary Classification: Active=1 vs Inactive=0)
   */
  class LogisticRegressionModel {
    constructor(learningRate = 0.05, iterations = 300, lambda = 0.01) {
      this.lr = learningRate;
      this.iterations = iterations;
      this.lambda = lambda;
      this.weights = [];
      this.bias = 0;
      this.featureNames = [];
    }

    sigmoid(z) {
      return 1 / (1 + Math.exp(-Math.max(-50, Math.min(50, z))));
    }

    fit(X, y, featureNames = []) {
      const N = X.length;
      const P = X[0].length;
      this.featureNames = featureNames;
      this.weights = new Array(P).fill(0);
      this.bias = 0;

      for (let iter = 0; iter < this.iterations; iter++) {
        const gradW = new Array(P).fill(0);
        let gradB = 0;

        for (let i = 0; i < N; i++) {
          let z = this.bias;
          for (let j = 0; j < P; j++) {
            z += this.weights[j] * X[i][j];
          }
          const pred = this.sigmoid(z);
          const err = pred - y[i];

          gradB += err;
          for (let j = 0; j < P; j++) {
            gradW[j] += err * X[i][j];
          }
        }

        // Apply gradient step with L2 regularization
        this.bias -= (this.lr * gradB) / N;
        for (let j = 0; j < P; j++) {
          this.weights[j] -= (this.lr * (gradW[j] + this.lambda * this.weights[j])) / N;
        }
      }
    }

    predictProba(X) {
      return X.map(row => {
        let z = this.bias;
        for (let j = 0; j < row.length; j++) {
          z += this.weights[j] * row[j];
        }
        return this.sigmoid(z);
      });
    }

    predict(X, threshold = 0.5) {
      return this.predictProba(X).map(p => (p >= threshold ? 1 : 0));
    }

    predictRow(row) {
      let z = this.bias;
      for (let j = 0; j < row.length; j++) {
        z += this.weights[j] * row[j];
      }
      const prob = this.sigmoid(z);
      return {
        probability: parseFloat(prob.toFixed(4)),
        prediction: prob >= 0.5 ? 1 : 0
      };
    }

    getFeatureImportance() {
      const total = this.weights.reduce((sum, w) => sum + Math.abs(w), 0) || 1.0;
      return this.weights.map((w, i) => ({
        feature: this.featureNames[i] || `Feature ${i + 1}`,
        weight: parseFloat(w.toFixed(4)),
        importance: parseFloat((Math.abs(w) / total).toFixed(4))
      })).sort((a, b) => b.importance - a.importance);
    }
  }

  /**
   * Calculate Regression Evaluation Metrics
   * @param {Array<number>} actual
   * @param {Array<number>} predicted
   * @returns {{ r2: number, rmse: number, mae: number, mse: number, residuals: Array<number> }}
   */
  function calculateRegressionMetrics(actual, predicted) {
    const N = actual.length;
    if (N === 0) return { r2: 0, rmse: 0, mae: 0, mse: 0, residuals: [] };

    let sumActual = 0;
    for (let i = 0; i < N; i++) {
      sumActual += actual[i];
    }
    const meanActual = sumActual / N;

    let ssTot = 0;
    let ssRes = 0;
    let sumAbsErr = 0;
    const residuals = [];

    for (let i = 0; i < N; i++) {
      const act = actual[i];
      const pred = predicted[i];
      const res = act - pred;
      residuals.push(res);

      ssRes += res * res;
      ssTot += (act - meanActual) * (act - meanActual);
      sumAbsErr += Math.abs(res);
    }

    const mse = ssRes / N;
    const rmse = Math.sqrt(mse);
    const mae = sumAbsErr / N;
    // R2 coefficient of determination (can be negative if model is worse than horizontal mean line)
    const r2 = ssTot === 0 ? 0 : Math.max(-1.0, 1 - ssRes / ssTot);

    return {
      r2: parseFloat(r2.toFixed(4)),
      rmse: parseFloat(rmse.toFixed(4)),
      mae: parseFloat(mae.toFixed(4)),
      mse: parseFloat(mse.toFixed(4)),
      residuals: residuals
    };
  }

  /**
   * Calculate Classification Evaluation Metrics
   * @param {Array<number>} actual
   * @param {Array<number>} predicted
   * @param {Array<number>} probabilities
   * @returns {{ accuracy: number, precision: number, recall: number, f1: number, confusionMatrix: Object }}
   */
  function calculateClassificationMetrics(actual, predicted, probabilities = []) {
    const N = actual.length;
    if (N === 0) {
      return {
        accuracy: 0,
        precision: 0,
        recall: 0,
        f1: 0,
        confusionMatrix: { tp: 0, tn: 0, fp: 0, fn: 0 }
      };
    }

    let tp = 0, tn = 0, fp = 0, fn = 0;
    for (let i = 0; i < N; i++) {
      const act = actual[i];
      const pred = predicted[i];
      if (act === 1 && pred === 1) tp++;
      else if (act === 0 && pred === 0) tn++;
      else if (act === 0 && pred === 1) fp++;
      else if (act === 1 && pred === 0) fn++;
    }

    const accuracy = (tp + tn) / N;
    const precision = tp + fp > 0 ? tp / (tp + fp) : 0;
    const recall = tp + fn > 0 ? tp / (tp + fn) : 0;
    const f1 = precision + recall > 0 ? (2 * precision * recall) / (precision + recall) : 0;

    return {
      accuracy: parseFloat(accuracy.toFixed(4)),
      precision: parseFloat(precision.toFixed(4)),
      recall: parseFloat(recall.toFixed(4)),
      f1: parseFloat(f1.toFixed(4)),
      confusionMatrix: { tp, tn, fp, fn }
    };
  }

  /**
   * K-Fold Cross Validation
   * @param {string} modelType - 'linear' or 'randomForest'
   * @param {Array} X
   * @param {Array} y
   * @param {number} k - folds (e.g. 5 or 10)
   * @param {Object} options
   */
  function kFoldCrossValidation(modelType, X, y, k = 5, options = {}) {
    const N = X.length;
    if (N < k * 2) {
      throw new Error(`Insufficient samples (${N}) for ${k}-fold cross validation.`);
    }

    const indices = Array.from({ length: N }, (_, i) => i);
    const prng = createPRNG(options.seed || 42);
    shuffleArray(indices, prng);

    const foldSize = Math.floor(N / k);
    const foldScores = [];

    for (let fold = 0; fold < k; fold++) {
      const valStart = fold * foldSize;
      const valEnd = fold === k - 1 ? N : valStart + foldSize;

      const valIndices = indices.slice(valStart, valEnd);
      const trainIndices = [...indices.slice(0, valStart), ...indices.slice(valEnd)];

      const XTrain = trainIndices.map(i => X[i]);
      const yTrain = trainIndices.map(i => y[i]);
      const XVal = valIndices.map(i => X[i]);
      const yVal = valIndices.map(i => y[i]);

      // Scale on fold training set
      const scaler = new StandardScaler();
      const XTrainScaled = scaler.fitTransform(XTrain);
      const XValScaled = scaler.transform(XVal);

      let model;
      if (modelType === 'randomForest') {
        model = new RandomForestRegressorModel(options.numTrees || 15, options.maxDepth || 5, 3, 0.7, (options.seed || 42) + fold);
      } else {
        model = new LinearRegressionModel(options.alpha || 1e-4);
      }

      model.fit(XTrainScaled, yTrain, options.featureNames);
      const yPred = model.predict(XValScaled);
      const metrics = calculateRegressionMetrics(yVal, yPred);
      foldScores.push(metrics.r2);
    }

    const meanScore = foldScores.reduce((a, b) => a + b, 0) / k;
    let variance = 0;
    foldScores.forEach(s => {
      variance += (s - meanScore) * (s - meanScore);
    });
    const stdScore = Math.sqrt(variance / k);

    return {
      k,
      foldScores: foldScores.map(s => parseFloat(s.toFixed(4))),
      meanScore: parseFloat(meanScore.toFixed(4)),
      stdScore: parseFloat(stdScore.toFixed(4))
    };
  }

  return {
    StandardScaler,
    trainTestSplit,
    LinearRegressionModel,
    RandomForestRegressorModel,
    LogisticRegressionModel,
    calculateRegressionMetrics,
    calculateClassificationMetrics,
    kFoldCrossValidation
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.ML = ML;
}
