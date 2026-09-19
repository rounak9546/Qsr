/**
 * qsar.js - QSAR Transformations, pIC50 Unit Handling, & Applicability Domain
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const QSAR = (function () {
  'use strict';

  // Unit conversion multiplier to Molar (M)
  const UNIT_MULTIPLIERS = {
    m: 1.0,
    mm: 1e-3,
    um: 1e-6,
    µm: 1e-6,
    nm: 1e-9,
    pm: 1e-12
  };

  /**
   * Convert an IC50 value with given unit to Molar concentration
   * @param {number} value
   * @param {string} unit (e.g. 'nM', 'uM', 'mM', 'M')
   * @returns {number|null} Molar concentration
   */
  function toMolar(value, unit) {
    if (!value || isNaN(value) || value <= 0) return null;
    const cleanUnit = (unit || 'nm').trim().toLowerCase();
    const mult = UNIT_MULTIPLIERS[cleanUnit];
    if (mult === undefined) return null;
    return value * mult;
  }

  /**
   * Calculate pIC50 from IC50 and its unit
   * Formula: pIC50 = -log10(IC50 in Molar)
   * Example: 10 nM = 1e-8 M -> pIC50 = -log10(1e-8) = 8.00
   * @param {number} value
   * @param {string} unit
   * @returns {{ pic50: number, molar: number, explanation: string } | null}
   */
  function calculatePIC50(value, unit) {
    const num = Number(value);
    if (isNaN(num) || num <= 0) return null;

    const cleanUnit = (unit || 'nM').trim();
    const molar = toMolar(num, cleanUnit);
    if (!molar || molar <= 0) return null;

    const pic50 = -Math.log10(molar);
    const roundedPic50 = parseFloat(pic50.toFixed(3));

    const explanation = `${num} ${cleanUnit} = ${molar.toExponential(2)} M → pIC₅₀ = -log₁₀(${molar.toExponential(2)}) = ${roundedPic50.toFixed(2)}`;

    return {
      pic50: roundedPic50,
      molar: molar,
      explanation: explanation
    };
  }

  /**
   * Convert pIC50 back to estimated IC50 in nM
   * Formula: IC50 (nM) = 10^(9 - pIC50)
   * @param {number} pic50
   * @returns {number}
   */
  function pic50ToNanomolar(pic50) {
    const p = Number(pic50);
    if (isNaN(p)) return 0;
    return Math.pow(10, 9 - p);
  }

  /**
   * Classify biological activity based on standard pharmaceutical criteria
   * Active: pIC50 >= 7.0 (IC50 <= 100 nM)
   * Intermediate: 6.0 <= pIC50 < 7.0 (100 nM < IC50 <= 1000 nM)
   * Inactive: pIC50 < 6.0 (IC50 > 1000 nM)
   * @param {number} pic50
   * @returns {{ label: string, color: string, badgeClass: string, binaryClass: number }}
   */
  function classifyActivity(pic50) {
    const p = Number(pic50);
    if (isNaN(p)) {
      return { label: 'Unknown', color: '#94a3b8', badgeClass: 'badge-unknown', binaryClass: 0 };
    }

    if (p >= 7.0) {
      return {
        label: 'Active (IC₅₀ ≤ 100 nM)',
        color: '#10b981', // emerald green
        badgeClass: 'badge-active',
        binaryClass: 1
      };
    } else if (p >= 6.0) {
      return {
        label: 'Intermediate (100 nM - 1 µM)',
        color: '#f59e0b', // amber
        badgeClass: 'badge-intermediate',
        binaryClass: 0
      };
    } else {
      return {
        label: 'Weak / Inactive (> 1 µM)',
        color: '#ef4444', // red
        badgeClass: 'badge-inactive',
        binaryClass: 0
      };
    }
  }

  /**
   * Compute Applicability Domain (AD) leverage and Euclidean distance in standardized feature space
   * @param {Array<Array<number>>} trainFeatures (standardized N x P)
   * @param {Array<number>} queryFeatures (standardized 1 x P)
   * @returns {{ insideDomain: boolean, distance: number, maxTrainDistance: number, message: string }}
   */
  function evaluateApplicabilityDomain(trainFeatures, queryFeatures) {
    if (!trainFeatures || trainFeatures.length === 0 || !queryFeatures) {
      return { insideDomain: true, distance: 0, maxTrainDistance: 0, message: 'No domain reference available.' };
    }

    const P = queryFeatures.length;
    // Calculate centroid of training features
    const centroid = new Array(P).fill(0);
    trainFeatures.forEach(row => {
      for (let j = 0; j < P; j++) {
        centroid[j] += row[j];
      }
    });
    for (let j = 0; j < P; j++) {
      centroid[j] /= trainFeatures.length;
    }

    // Compute training set distances to centroid
    let maxDist = 0;
    let distSum = 0;
    const trainDistances = trainFeatures.map(row => {
      let sumSq = 0;
      for (let j = 0; j < P; j++) {
        const d = row[j] - centroid[j];
        sumSq += d * d;
      }
      const dist = Math.sqrt(sumSq);
      if (dist > maxDist) maxDist = dist;
      distSum += dist;
      return dist;
    });

    const meanDist = distSum / trainFeatures.length;
    // Standard deviation of distances
    let variance = 0;
    trainDistances.forEach(d => {
      variance += (d - meanDist) * (d - meanDist);
    });
    const sdDist = Math.sqrt(variance / trainFeatures.length);
    const threshold = meanDist + 2.5 * sdDist; // 99% coverage boundary

    // Query distance to centroid
    let qSumSq = 0;
    for (let j = 0; j < P; j++) {
      const d = queryFeatures[j] - centroid[j];
      qSumSq += d * d;
    }
    const queryDist = Math.sqrt(qSumSq);

    const inside = queryDist <= threshold;
    const message = inside
      ? `Within Applicability Domain (Distance: ${queryDist.toFixed(2)} ≤ Threshold: ${threshold.toFixed(2)})`
      : `Extrapolation Warning: Candidate lies outside applicability domain (Distance: ${queryDist.toFixed(2)} > Threshold: ${threshold.toFixed(2)}). Prediction may have higher uncertainty.`;

    return {
      insideDomain: inside,
      distance: parseFloat(queryDist.toFixed(3)),
      threshold: parseFloat(threshold.toFixed(3)),
      message: message
    };
  }

  return {
    UNIT_MULTIPLIERS,
    toMolar,
    calculatePIC50,
    pic50ToNanomolar,
    classifyActivity,
    evaluateApplicabilityDomain
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.QSAR = QSAR;
}
