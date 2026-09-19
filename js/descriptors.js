/**
 * descriptors.js - Molecular Descriptors Calculation & Chemical Structure Renderer
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const Descriptors = (function () {
  'use strict';

  // Standard monoisotopic/average atomic weights (g/mol)
  const ATOMIC_WEIGHTS = {
    H: 1.008,
    C: 12.011,
    N: 14.007,
    O: 15.999,
    F: 18.998,
    P: 30.974,
    S: 32.065,
    Cl: 35.453,
    Br: 79.904,
    I: 126.904,
    B: 10.811
  };

  /**
   * Parse a SMILES string into a graph of atoms and bonds
   * @param {string} smiles
   * @returns {{ atoms: Array, bonds: Array, rings: number, aromaticRings: number }}
   */
  function parseSMILES(smiles) {
    if (!smiles || typeof smiles !== 'string') {
      return { atoms: [], bonds: [], rings: 0, aromaticRings: 0 };
    }

    const atoms = [];
    const bonds = [];
    const ringOpenings = {};
    const branchStack = [];
    let currentAtomIndex = -1;
    let pendingBondOrder = 1; // default single bond

    let i = 0;
    while (i < smiles.length) {
      const char = smiles[i];

      if (char === '(') {
        branchStack.push(currentAtomIndex);
        i++;
        continue;
      }

      if (char === ')') {
        if (branchStack.length > 0) {
          currentAtomIndex = branchStack.pop();
        }
        i++;
        continue;
      }

      // Bond specifications
      if (char === '-') {
        pendingBondOrder = 1;
        i++;
        continue;
      }
      if (char === '=') {
        pendingBondOrder = 2;
        i++;
        continue;
      }
      if (char === '#') {
        pendingBondOrder = 3;
        i++;
        continue;
      }
      if (char === ':') {
        pendingBondOrder = 1.5; // aromatic bond
        i++;
        continue;
      }
      if (char === '/' || char === '\\' || char === '.') {
        // stereochemistry or disconnect
        i++;
        continue;
      }

      // Explicit atom inside brackets [ ... ]
      if (char === '[') {
        const closeIdx = smiles.indexOf(']', i);
        if (closeIdx !== -1) {
          const content = smiles.substring(i + 1, closeIdx);
          const atom = parseBracketAtom(content, atoms.length);
          atoms.push(atom);
          const newIdx = atoms.length - 1;

          if (currentAtomIndex !== -1) {
            bonds.push({
              source: currentAtomIndex,
              target: newIdx,
              order: pendingBondOrder,
              isAromatic: atom.isAromatic && (atoms[currentAtomIndex].isAromatic || false)
            });
            pendingBondOrder = 1;
          }
          currentAtomIndex = newIdx;
          i = closeIdx + 1;
          continue;
        }
      }

      // Ring closure digit: 1-9 or %digit
      if ((char >= '1' && char <= '9') || char === '%') {
        let ringNum;
        if (char === '%') {
          ringNum = smiles.substring(i + 1, i + 3);
          i += 3;
        } else {
          ringNum = char;
          i++;
        }

        if (ringOpenings[ringNum] !== undefined) {
          const partnerIdx = ringOpenings[ringNum];
          bonds.push({
            source: partnerIdx,
            target: currentAtomIndex,
            order: pendingBondOrder,
            isRing: true
          });
          delete ringOpenings[ringNum];
          pendingBondOrder = 1;
        } else {
          ringOpenings[ringNum] = currentAtomIndex;
        }
        continue;
      }

      // Standard organic subset atoms: Cl, Br, C, c, N, n, O, o, S, s, P, p, F, I, B
      let symbol = null;
      let isAromatic = false;

      if (i + 1 < smiles.length && (smiles.substring(i, i + 2) === 'Cl' || smiles.substring(i, i + 2) === 'Br')) {
        symbol = smiles.substring(i, i + 2);
        i += 2;
      } else if ('BCNOPSFIbcnops'.includes(char)) {
        symbol = char;
        if (char >= 'a' && char <= 'z') {
          isAromatic = true;
          symbol = char.toUpperCase();
        }
        i++;
      } else {
        // Unknown symbol or skip
        i++;
        continue;
      }

      if (symbol) {
        const atom = {
          index: atoms.length,
          symbol: symbol,
          isAromatic: isAromatic,
          charge: 0,
          explicitH: 0
        };
        atoms.push(atom);
        const newIdx = atoms.length - 1;

        if (currentAtomIndex !== -1) {
          bonds.push({
            source: currentAtomIndex,
            target: newIdx,
            order: (isAromatic && atoms[currentAtomIndex].isAromatic) ? 1.5 : pendingBondOrder,
            isAromatic: isAromatic && atoms[currentAtomIndex].isAromatic
          });
          pendingBondOrder = 1;
        }
        currentAtomIndex = newIdx;
      }
    }

    // Estimate rings from Euler-Poincaré formula for planar graph: R = Bonds - Atoms + ConnectedComponents
    const bondCount = bonds.length;
    const atomCount = atoms.length;
    const estimatedRings = Math.max(0, bondCount - atomCount + 1);

    // Count aromatic rings: groups of ~5-6 connected aromatic atoms
    const aromaticAtoms = atoms.filter(a => a.isAromatic).length;
    const aromaticRings = Math.floor(aromaticAtoms / 5.5);

    return {
      atoms,
      bonds,
      rings: estimatedRings,
      aromaticRings: Math.min(estimatedRings, aromaticRings)
    };
  }

  /**
   * Helper to parse bracketed atom: e.g. [NH+], [O-], [nH]
   */
  function parseBracketAtom(content, index) {
    let symbol = 'C';
    let isAromatic = false;
    let explicitH = 0;
    let charge = 0;

    // Detect charge: + or -
    if (content.includes('+')) {
      const match = content.match(/\+(\d*)/);
      charge = match && match[1] ? parseInt(match[1], 10) : 1;
    } else if (content.includes('-')) {
      const match = content.match(/\-(\d*)/);
      charge = match && match[1] ? -parseInt(match[1], 10) : -1;
    }

    // Detect explicit H: H or H2, H3
    if (content.includes('H')) {
      const match = content.match(/H(\d*)/);
      explicitH = match && match[1] ? parseInt(match[1], 10) : 1;
    }

    // Extract elemental symbol
    const cleanContent = content.replace(/[\+\-\d@H]/g, '');
    if (cleanContent.length > 0) {
      if (cleanContent.length >= 2 && ['Cl', 'Br', 'Si', 'Na', 'Li'].includes(cleanContent.substring(0, 2))) {
        symbol = cleanContent.substring(0, 2);
      } else {
        const firstChar = cleanContent[0];
        if (firstChar >= 'a' && firstChar <= 'z') {
          isAromatic = true;
          symbol = firstChar.toUpperCase();
        } else {
          symbol = firstChar;
        }
      }
    }

    return {
      index,
      symbol,
      isAromatic,
      charge,
      explicitH
    };
  }

  /**
   * Calculate 2D coordinates for chemical visualization using spring-electrical force layout
   * @param {Array} atoms
   * @param {Array} bonds
   * @param {number} width
   * @param {number} height
   */
  function layoutMolecule(atoms, bonds, width, height) {
    const N = atoms.length;
    if (N === 0) return [];

    const coords = [];
    const radius = Math.min(width, height) * 0.35;
    const cx = width / 2;
    const cy = height / 2;

    // Initial positioning in circular arrangement
    for (let i = 0; i < N; i++) {
      const angle = (2 * Math.PI * i) / N;
      coords.push({
        x: cx + radius * Math.cos(angle) + (Math.random() - 0.5) * 10,
        y: cy + radius * Math.sin(angle) + (Math.random() - 0.5) * 10,
        vx: 0,
        vy: 0
      });
    }

    // Force-directed relaxation (spring-electrical embedding)
    const iterations = Math.min(80, 20 + N * 2);
    const kTarget = 38; // desired bond length in pixels
    const kRepel = 700;

    for (let iter = 0; iter < iterations; iter++) {
      const temp = 1.0 - iter / iterations;

      // Coulomb-like repulsive forces between all atom pairs
      for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
          const dx = coords[i].x - coords[j].x;
          const dy = coords[i].y - coords[j].y;
          const distSq = dx * dx + dy * dy + 0.01;
          const dist = Math.sqrt(distSq);

          const force = (kRepel * kRepel) / distSq;
          const fx = (dx / dist) * force;
          const fy = (dy / dist) * force;

          coords[i].vx += fx;
          coords[i].vy += fy;
          coords[j].vx -= fx;
          coords[j].vy -= fy;
        }
      }

      // Hooke-like spring attractive forces along bonds
      for (let b = 0; b < bonds.length; b++) {
        const bond = bonds[b];
        const i = bond.source;
        const j = bond.target;
        if (i >= N || j >= N) continue;

        const dx = coords[j].x - coords[i].x;
        const dy = coords[j].y - coords[i].y;
        const dist = Math.sqrt(dx * dx + dy * dy) + 0.001;
        const disp = dist - kTarget;

        const force = 0.08 * disp;
        const fx = (dx / dist) * force;
        const fy = (dy / dist) * force;

        coords[i].vx += fx;
        coords[i].vy += fy;
        coords[j].vx -= fx;
        coords[j].vy -= fy;
      }

      // Update positions with damping
      for (let i = 0; i < N; i++) {
        coords[i].x += coords[i].vx * 0.25 * temp;
        coords[i].y += coords[i].vy * 0.25 * temp;
        coords[i].vx *= 0.6;
        coords[i].vy *= 0.6;
      }
    }

    // Normalize coordinates within canvas viewport with comfortable margin
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (let i = 0; i < N; i++) {
      minX = Math.min(minX, coords[i].x);
      maxX = Math.max(maxX, coords[i].x);
      minY = Math.min(minY, coords[i].y);
      maxY = Math.max(maxY, coords[i].y);
    }

    const padding = 35;
    const molW = Math.max(10, maxX - minX);
    const molH = Math.max(10, maxY - minY);
    const scale = Math.min((width - 2 * padding) / molW, (height - 2 * padding) / molH, 1.4);

    const midX = (minX + maxX) / 2;
    const midY = (minY + maxY) / 2;

    for (let i = 0; i < N; i++) {
      coords[i].x = cx + (coords[i].x - midX) * scale;
      coords[i].y = cy + (coords[i].y - midY) * scale;
    }

    return coords;
  }

  /**
   * Render molecular graph on HTML5 Canvas
   * @param {HTMLCanvasElement} canvas
   * @param {string} smiles
   * @param {Object} [options]
   */
  function renderMolecule(canvas, smiles, options = {}) {
    if (!canvas || !canvas.getContext) return false;
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;

    // Clear canvas
    ctx.clearRect(0, 0, width, height);

    // Background styling
    const bgColor = options.bgColor || '#0f172a';
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, width, height);

    const mol = parseSMILES(smiles);
    if (mol.atoms.length === 0) {
      ctx.fillStyle = '#94a3b8';
      ctx.font = '14px Inter, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('Structure unavailable or empty SMILES', width / 2, height / 2);
      return false;
    }

    const coords = layoutMolecule(mol.atoms, mol.bonds, width, height);

    // Color mapping for heteroatoms
    const atomColors = {
      C: '#cbd5e1', // Slate light
      N: '#38bdf8', // Cyan/sky blue
      O: '#f87171', // Coral red
      S: '#fbbf24', // Amber
      F: '#34d399', // Emerald
      Cl: '#4ade80', // Green
      Br: '#a855f7', // Purple
      I: '#e879f9', // Magenta
      P: '#fb923c'  // Orange
    };

    // Draw Bonds
    ctx.lineWidth = 2.2;
    ctx.lineCap = 'round';

    mol.bonds.forEach(bond => {
      const p1 = coords[bond.source];
      const p2 = coords[bond.target];
      if (!p1 || !p2) return;

      const dx = p2.x - p1.x;
      const dy = p2.y - p1.y;
      const dist = Math.sqrt(dx * dx + dy * dy);
      const nx = -dy / dist;
      const ny = dx / dist;

      if (bond.order === 2) {
        // Double bond
        const offset = 2.4;
        ctx.strokeStyle = '#64748b';
        ctx.beginPath();
        ctx.moveTo(p1.x + nx * offset, p1.y + ny * offset);
        ctx.lineTo(p2.x + nx * offset, p2.y + ny * offset);
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(p1.x - nx * offset, p1.y - ny * offset);
        ctx.lineTo(p2.x - nx * offset, p2.y - ny * offset);
        ctx.stroke();
      } else if (bond.order === 3) {
        // Triple bond
        const offset = 3.6;
        ctx.strokeStyle = '#64748b';
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(p1.x + nx * offset, p1.y + ny * offset);
        ctx.lineTo(p2.x + nx * offset, p2.y + ny * offset);
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(p1.x - nx * offset, p1.y - ny * offset);
        ctx.lineTo(p2.x - nx * offset, p2.y - ny * offset);
        ctx.stroke();
      } else if (bond.isAromatic) {
        // Aromatic bond (dashed style)
        ctx.strokeStyle = '#38bdf8';
        ctx.setLineDash([4, 3]);
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.stroke();
        ctx.setLineDash([]);
      } else {
        // Single bond
        ctx.strokeStyle = '#64748b';
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.stroke();
      }
    });

    // Draw Atoms
    mol.atoms.forEach(atom => {
      const pt = coords[atom.index];
      if (!pt) return;

      const isCarbon = atom.symbol === 'C' && !atom.isAromatic;
      const color = atomColors[atom.symbol] || '#94a3b8';

      if (!isCarbon || atom.charge !== 0) {
        // Draw heteroatom circle background
        ctx.fillStyle = bgColor;
        ctx.beginPath();
        ctx.arc(pt.x, pt.y, 9, 0, 2 * Math.PI);
        ctx.fill();

        // Draw symbol
        ctx.fillStyle = color;
        ctx.font = 'bold 12px "JetBrains Mono", Consolas, monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';

        let label = atom.symbol;
        if (atom.charge === 1) label += '+';
        if (atom.charge === -1) label += '⁻';
        ctx.fillText(label, pt.x, pt.y);
      } else {
        // Small node dot for carbon vertices in dark mode
        ctx.fillStyle = '#475569';
        ctx.beginPath();
        ctx.arc(pt.x, pt.y, 2.2, 0, 2 * Math.PI);
        ctx.fill();
      }
    });

    return true;
  }

  /**
   * Estimate molecular descriptors directly from SMILES chemical graph
   * @param {string} smiles
   * @returns {Object} Descriptors: { molecular_weight, logP, hbd, hba, rotatable_bonds, tpsa, ring_count, aromatic_ring_count, fraction_sp3, lipinski_violations }
   */
  function calculateDescriptorsFromSMILES(smiles) {
    const mol = parseSMILES(smiles);
    const atoms = mol.atoms;
    const bonds = mol.bonds;

    let mw = 0;
    let hbd = 0;
    let hba = 0;
    let tpsa = 0;
    let sp3Carbons = 0;
    let totalCarbons = 0;
    let rotatableBonds = 0;

    // Hydrogen counts estimation based on standard valences
    const valenceMap = { C: 4, N: 3, O: 2, S: 2, P: 3, F: 1, Cl: 1, Br: 1, I: 1 };

    // TPSA fragment contributions (Ertl et al., 2000)
    // -OH: ~20.23, -NH2: ~26.02, -NH-: ~12.03, =N-: ~12.89, -O-: ~9.23, =O: ~17.07
    atoms.forEach(atom => {
      const sym = atom.symbol;
      const weight = ATOMIC_WEIGHTS[sym] || 12.0;

      // Connected bond order sum
      let bondOrderSum = 0;
      let connectedNeighbors = 0;
      bonds.forEach(b => {
        if (b.source === atom.index || b.target === atom.index) {
          bondOrderSum += b.order;
          connectedNeighbors++;
        }
      });

      const maxValence = valenceMap[sym] || 4;
      const implicitH = Math.max(0, maxValence - Math.round(bondOrderSum) + atom.charge);
      const totalH = atom.explicitH > 0 ? atom.explicitH : implicitH;

      mw += weight + totalH * ATOMIC_WEIGHTS.H;

      if (sym === 'C') {
        totalCarbons++;
        if (!atom.isAromatic && bondOrderSum <= 4 && connectedNeighbors + totalH === 4) {
          sp3Carbons++;
        }
      }

      // Hydrogen Bond Donors (OH, NH, NH2)
      if ((sym === 'N' || sym === 'O') && totalH > 0) {
        hbd += totalH;
      }

      // Hydrogen Bond Acceptors (N, O with lone pairs)
      if (sym === 'O') {
        hba += 1;
        // TPSA contribution for Oxygen
        if (totalH > 0) tpsa += 20.23; // -OH
        else if (bondOrderSum >= 2 && connectedNeighbors === 1) tpsa += 17.07; // =O
        else tpsa += 9.23; // -O-
      } else if (sym === 'N') {
        hba += 1;
        // TPSA contribution for Nitrogen
        if (totalH === 2) tpsa += 26.02; // -NH2
        else if (totalH === 1) tpsa += 12.03; // -NH-
        else if (bondOrderSum >= 3 && connectedNeighbors === 1) tpsa += 23.79; // #N
        else tpsa += 12.89; // =N- or >N-
      } else if (sym === 'S') {
        // Polar sulfur
        if (connectedNeighbors >= 3) tpsa += 25.0; // sulfone/sulfoxide
      }
    });

    // Rotatable bonds: single non-ring bonds between heavy non-terminal atoms
    bonds.forEach(b => {
      if (b.order === 1 && !b.isRing && !b.isAromatic) {
        const a1 = atoms[b.source];
        const a2 = atoms[b.target];
        if (a1 && a2 && a1.symbol !== 'H' && a2.symbol !== 'H') {
          // Check terminal atom status
          const n1 = bonds.filter(x => x.source === b.source || x.target === b.source).length;
          const n2 = bonds.filter(x => x.source === b.target || x.target === b.target).length;
          if (n1 > 1 && n2 > 1) {
            rotatableBonds++;
          }
        }
      }
    });

    // Wildman-Crippen LogP estimation approximation based on lipophilic and hydrophilic contributions
    // Hydrophobic core: C (+0.25), aromatic C (+0.35), Cl (+0.6), Br (+0.8), F (+0.15)
    // Hydrophilic: N (-0.3), O (-0.4), OH (-0.7), Polar rings/amide (-0.5)
    let estimatedLogP = 0.5;
    atoms.forEach(a => {
      if (a.symbol === 'C') estimatedLogP += a.isAromatic ? 0.35 : 0.22;
      else if (a.symbol === 'N') estimatedLogP -= 0.32;
      else if (a.symbol === 'O') estimatedLogP -= 0.45;
      else if (a.symbol === 'F') estimatedLogP += 0.15;
      else if (a.symbol === 'Cl') estimatedLogP += 0.65;
      else if (a.symbol === 'Br') estimatedLogP += 0.85;
      else if (a.symbol === 'I') estimatedLogP += 1.10;
      else if (a.symbol === 'S') estimatedLogP += 0.30;
    });

    const fractionSp3 = totalCarbons > 0 ? sp3Carbons / totalCarbons : 0;

    // Evaluate Lipinski's Rule of 5 Violations
    let lipinskiViolations = 0;
    const lipinskiDetails = [];

    if (mw > 500) {
      lipinskiViolations++;
      lipinskiDetails.push(`MW > 500 (${mw.toFixed(1)} Da)`);
    }
    if (estimatedLogP > 5.0) {
      lipinskiViolations++;
      lipinskiDetails.push(`LogP > 5.0 (${estimatedLogP.toFixed(2)})`);
    }
    if (hbd > 5) {
      lipinskiViolations++;
      lipinskiDetails.push(`HBD > 5 (${hbd})`);
    }
    if (hba > 10) {
      lipinskiViolations++;
      lipinskiDetails.push(`HBA > 10 (${hba})`);
    }

    // Veber Rule (RotBonds <= 10, TPSA <= 140 Å²)
    const veberCompliant = rotatableBonds <= 10 && tpsa <= 140;

    return {
      molecular_weight: parseFloat(mw.toFixed(2)),
      logP: parseFloat(estimatedLogP.toFixed(2)),
      hbd: Math.round(hbd),
      hba: Math.round(hba),
      rotatable_bonds: Math.round(rotatableBonds),
      tpsa: parseFloat(tpsa.toFixed(2)),
      ring_count: mol.rings,
      aromatic_ring_count: mol.aromaticRings,
      fraction_sp3: parseFloat(fractionSp3.toFixed(2)),
      lipinski_violations: lipinskiViolations,
      lipinski_details: lipinskiDetails,
      veber_compliant: veberCompliant,
      provenance: 'Computed in-browser (Deterministic Chemoinformatics Engine)'
    };
  }

  return {
    parseSMILES,
    renderMolecule,
    calculateDescriptorsFromSMILES,
    ATOMIC_WEIGHTS
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.Descriptors = Descriptors;
}
