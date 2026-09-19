/**
 * validation.js - Data, SMILES, and Input Validation Module
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const Validator = (function () {
  'use strict';

  /**
   * Sanitize text string for safe HTML rendering
   * @param {string} str
   * @returns {string}
   */
  function sanitize(str) {
    if (str === null || str === undefined) return '';
    return String(str)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#039;');
  }

  /**
   * Validate a SMILES string using structural heuristics
   * Checks for balanced brackets, valid chemical symbols, valid ring closures, and non-empty input
   * @param {string} smiles
   * @returns {{ valid: boolean, error?: string, normalized?: string }}
   */
  function validateSMILES(smiles) {
    if (!smiles || typeof smiles !== 'string') {
      return { valid: false, error: 'SMILES string is empty or invalid.' };
    }

    const trimmed = smiles.trim();
    if (trimmed.length === 0) {
      return { valid: false, error: 'SMILES string cannot be empty.' };
    }

    // Check for invalid characters (standard SMILES permits organic subset, brackets, bonds, rings, charges)
    // Allowed characters: C, N, O, S, P, F, Cl, Br, I, B, H, c, n, o, s, p, brackets [], parens (), bonds - = # : / \ . @ + -, numbers 0-9, %
    const validCharsRegex = /^[A-Za-z0-9\(\)\[\]\=\#\:\-\+\/\\\.\@\%\*\$]+$/;
    if (!validCharsRegex.test(trimmed)) {
      return { valid: false, error: 'SMILES contains invalid chemical characters or unsupported symbols.' };
    }

    // Check balanced parentheses (branching)
    let parenCount = 0;
    for (let i = 0; i < trimmed.length; i++) {
      if (trimmed[i] === '(') parenCount++;
      if (trimmed[i] === ')') parenCount--;
      if (parenCount < 0) {
        return { valid: false, error: 'Unbalanced parentheses in SMILES branch structure.' };
      }
    }
    if (parenCount !== 0) {
      return { valid: false, error: 'Unbalanced parentheses in SMILES branch structure.' };
    }

    // Check balanced square brackets (explicit atoms)
    let bracketCount = 0;
    for (let i = 0; i < trimmed.length; i++) {
      if (trimmed[i] === '[') bracketCount++;
      if (trimmed[i] === ']') bracketCount--;
      if (bracketCount < 0) {
        return { valid: false, error: 'Unbalanced square brackets in SMILES atom definition.' };
      }
    }
    if (bracketCount !== 0) {
      return { valid: false, error: 'Unbalanced square brackets in SMILES atom definition.' };
    }

    // Check ring closure pairing (digits 1-9 should occur in pairs, unless multi-cycle)
    const digitCounts = {};
    let inBracket = false;
    for (let i = 0; i < trimmed.length; i++) {
      const ch = trimmed[i];
      if (ch === '[') inBracket = true;
      else if (ch === ']') inBracket = false;
      else if (!inBracket && ch >= '1' && ch <= '9') {
        digitCounts[ch] = (digitCounts[ch] || 0) + 1;
      }
    }

    for (const [digit, count] of Object.entries(digitCounts)) {
      if (count % 2 !== 0) {
        return {
          valid: false,
          error: `Unmatched ring closure number '${digit}' found (${count} occurrences). SMILES ring closures must occur in pairs.`
        };
      }
    }

    return { valid: true, normalized: trimmed };
  }

  /**
   * Validate CSV columns for required fields
   * @param {Array<string>} headers
   * @returns {{ valid: boolean, missing: Array<string> }}
   */
  function validateCSVHeaders(headers) {
    if (!Array.isArray(headers) || headers.length === 0) {
      return { valid: false, missing: ['headers'] };
    }

    const normalizedHeaders = headers.map(h => (h ? h.trim().toLowerCase().replace(/[\s_-]+/g, '') : ''));

    // Minimum essential columns for QSAR: SMILES, Activity Value or pIC50
    const hasSmiles = normalizedHeaders.some(h => h.includes('smiles'));
    const hasActivity = normalizedHeaders.some(h =>
      h.includes('activity') || h.includes('ic50') || h.includes('pic50') || h.includes('value')
    );

    const missing = [];
    if (!hasSmiles) missing.push('SMILES column');
    if (!hasActivity) missing.push('Activity value or IC50 / pIC50 column');

    return {
      valid: missing.length === 0,
      missing: missing
    };
  }

  /**
   * Check if a value is a valid, finite positive number
   * @param {*} val
   * @returns {boolean}
   */
  function isPositiveNumber(val) {
    const num = Number(val);
    return !isNaN(num) && isFinite(num) && num > 0;
  }

  /**
   * Validate numeric range for descriptors
   * @param {Object} desc
   * @returns {Array<string>} list of warnings/errors
   */
  function validateDescriptorRanges(desc) {
    const issues = [];
    if (desc.molecular_weight !== undefined && (desc.molecular_weight <= 0 || desc.molecular_weight > 2500)) {
      issues.push(`Molecular weight (${desc.molecular_weight}) is outside realistic small molecule range (0-2500 Da).`);
    }
    if (desc.logP !== undefined && (desc.logP < -10 || desc.logP > 15)) {
      issues.push(`LogP (${desc.logP}) is outside typical range (-10 to 15).`);
    }
    if (desc.tpsa !== undefined && (desc.tpsa < 0 || desc.tpsa > 800)) {
      issues.push(`TPSA (${desc.tpsa}) is outside typical range (0-800 Å²).`);
    }
    return issues;
  }

  return {
    sanitize,
    validateSMILES,
    validateCSVHeaders,
    isPositiveNumber,
    validateDescriptorRanges
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.Validator = Validator;
}
