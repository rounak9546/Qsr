/**
 * ui.js - User Interface Controller, Navigation, Modals, Tables, and Toasts
 * Academic Bioinformatics / Computational Drug Discovery Platform
 */

const UI = (function () {
  'use strict';

  // Navigation state
  let activeSectionId = 'overview';

  // Table pagination state for Dataset Explorer
  let datasetTableState = {
    page: 1,
    pageSize: 10,
    searchQuery: '',
    activityFilter: 'all',
    sortColumn: 'pic50',
    sortDirection: 'desc'
  };

  // Table pagination state for Candidate Prioritization
  let candidateTableState = {
    page: 1,
    pageSize: 10,
    searchQuery: '',
    sortColumn: 'predicted_pic50',
    sortDirection: 'desc'
  };

  /**
   * Show Toast Notification in top-right corner
   * @param {string} message
   * @param {string} type - 'success' | 'warning' | 'error' | 'info'
   * @param {number} duration - ms
   */
  function showToast(message, type = 'info', duration = 4000) {
    const container = document.getElementById('toast-container');
    if (!container) return;

    const toast = document.createElement('div');
    toast.className = `toast toast-${type}`;
    
    // Icon based on type
    let icon = 'ℹ️';
    if (type === 'success') icon = '✓';
    if (type === 'warning') icon = '⚠️';
    if (type === 'error') icon = '✕';

    toast.innerHTML = `
      <div class="toast-icon">${icon}</div>
      <div class="toast-content">${Validator.sanitize(message)}</div>
      <button class="toast-close" aria-label="Close notification">&times;</button>
    `;

    container.appendChild(toast);

    // Auto-remove
    const timer = setTimeout(() => {
      removeToast(toast);
    }, duration);

    toast.querySelector('.toast-close').addEventListener('click', () => {
      clearTimeout(timer);
      removeToast(toast);
    });
  }

  function removeToast(toast) {
    toast.classList.add('toast-fadeout');
    setTimeout(() => {
      if (toast.parentNode) {
        toast.parentNode.removeChild(toast);
      }
    }, 300);
  }

  /**
   * Switch Active Section in SPA Navigation
   * @param {string} sectionId
   */
  function navigateTo(sectionId) {
    const target = document.getElementById(sectionId);
    if (!target) {
      console.warn(`Section #${sectionId} not found.`);
      return;
    }

    // Hide all sections
    const sections = document.querySelectorAll('.app-section');
    sections.forEach(sec => {
      sec.classList.remove('active');
    });

    // Show target section
    target.classList.add('active');
    activeSectionId = sectionId;

    // Update active state on navigation buttons/links
    const navLinks = document.querySelectorAll('.nav-item');
    navLinks.forEach(link => {
      if (link.getAttribute('data-target') === sectionId) {
        link.classList.add('active');
      } else {
        link.classList.remove('active');
      }
    });

    // Close mobile drawer if open
    const sidebar = document.getElementById('app-sidebar');
    if (sidebar && sidebar.classList.contains('mobile-open')) {
      sidebar.classList.remove('mobile-open');
    }

    // Scroll to top of main content container smoothly
    const mainContent = document.getElementById('main-content');
    if (mainContent) {
      mainContent.scrollTo({ top: 0, behavior: 'smooth' });
    }

    // Update URL hash without reload
    if (window.history.replaceState) {
      window.history.replaceState(null, null, `#${sectionId}`);
    }

    // Trigger chart resize / redraw if navigating to analysis or results
    window.dispatchEvent(new Event('resize'));
  }

  /**
   * Render Dataset Explorer Table with search, sorting, filtering, and pagination
   * @param {Array<Object>} dataset
   */
  function renderDatasetTable(dataset) {
    const tableBody = document.getElementById('dataset-table-body');
    const paginationInfo = document.getElementById('dataset-pagination-info');
    const prevBtn = document.getElementById('dataset-prev-page');
    const nextBtn = document.getElementById('dataset-next-page');
    const countBadge = document.getElementById('dataset-count-badge');

    if (!tableBody) return;

    if (!dataset || dataset.length === 0) {
      tableBody.innerHTML = `
        <tr>
          <td colspan="10" class="text-center py-8 text-muted">
            No compounds in dataset. Upload a CSV or load the benchmark dataset.
          </td>
        </tr>
      `;
      if (countBadge) countBadge.textContent = '0 compounds';
      return;
    }

    // Filter by search query
    let filtered = dataset.filter(item => {
      const q = datasetTableState.searchQuery.toLowerCase();
      if (!q) return true;
      return (
        (item.compound_id && item.compound_id.toLowerCase().includes(q)) ||
        (item.compound_name && item.compound_name.toLowerCase().includes(q)) ||
        (item.smiles && item.smiles.toLowerCase().includes(q))
      );
    });

    // Filter by activity class
    if (datasetTableState.activityFilter !== 'all') {
      filtered = filtered.filter(item => {
        if (datasetTableState.activityFilter === 'active') return item.pic50 >= 7.0;
        if (datasetTableState.activityFilter === 'intermediate') return item.pic50 >= 6.0 && item.pic50 < 7.0;
        if (datasetTableState.activityFilter === 'inactive') return item.pic50 < 6.0;
        return true;
      });
    }

    // Sort
    const col = datasetTableState.sortColumn;
    const dir = datasetTableState.sortDirection === 'asc' ? 1 : -1;
    filtered.sort((a, b) => {
      let vA = a[col];
      let vB = b[col];
      if (typeof vA === 'string') {
        return vA.localeCompare(vB) * dir;
      }
      return ((Number(vA) || 0) - (Number(vB) || 0)) * dir;
    });

    if (countBadge) {
      countBadge.textContent = `${filtered.length} of ${dataset.length} compounds`;
    }

    // Pagination slice
    const totalPages = Math.ceil(filtered.length / datasetTableState.pageSize) || 1;
    if (datasetTableState.page > totalPages) datasetTableState.page = totalPages;
    if (datasetTableState.page < 1) datasetTableState.page = 1;

    const startIdx = (datasetTableState.page - 1) * datasetTableState.pageSize;
    const pageItems = filtered.slice(startIdx, startIdx + datasetTableState.pageSize);

    if (pageItems.length === 0) {
      tableBody.innerHTML = `
        <tr>
          <td colspan="10" class="text-center py-6 text-muted">
            No matching compounds found for search query "${Validator.sanitize(datasetTableState.searchQuery)}".
          </td>
        </tr>
      `;
    } else {
      tableBody.innerHTML = pageItems
        .map(item => {
          let badgeColor = '#ef4444';
          let badgeText = 'Inactive';
          if (item.pic50 >= 7.0) {
            badgeColor = '#10b981';
            badgeText = 'Active';
          } else if (item.pic50 >= 6.0) {
            badgeColor = '#f59e0b';
            badgeText = 'Intermediate';
          }

          return `
          <tr data-compound-id="${Validator.sanitize(item.compound_id)}">
            <td><strong class="text-cyan">${Validator.sanitize(item.compound_id)}</strong></td>
            <td>${Validator.sanitize(item.compound_name)}</td>
            <td class="font-mono text-xs max-w-xs truncate" title="${Validator.sanitize(item.smiles)}">${Validator.sanitize(item.smiles)}</td>
            <td><strong>${item.pic50.toFixed(2)}</strong></td>
            <td>${item.activity_value.toFixed(1)} ${item.activity_unit || 'nM'}</td>
            <td><span class="badge" style="background-color: ${badgeColor}22; color: ${badgeColor}; border: 1px solid ${badgeColor}44;">${badgeText}</span></td>
            <td>${item.molecular_weight.toFixed(1)}</td>
            <td>${item.logP.toFixed(2)}</td>
            <td>${item.tpsa.toFixed(1)}</td>
            <td>
              <button class="btn-inspect btn-sm" onclick="UI.showCompoundModal('${Validator.sanitize(item.compound_id)}')">
                Inspect
              </button>
            </td>
          </tr>
        `;
        })
        .join('');
    }

    if (paginationInfo) {
      paginationInfo.textContent = `Page ${datasetTableState.page} of ${totalPages} (${filtered.length} total)`;
    }
    if (prevBtn) prevBtn.disabled = datasetTableState.page <= 1;
    if (nextBtn) nextBtn.disabled = datasetTableState.page >= totalPages;
  }

  /**
   * Render Candidate Ranking Table
   * @param {Array<Object>} candidates
   */
  function renderCandidateTable(candidates) {
    const tableBody = document.getElementById('candidate-table-body');
    const paginationInfo = document.getElementById('candidate-pagination-info');
    const prevBtn = document.getElementById('candidate-prev-page');
    const nextBtn = document.getElementById('candidate-next-page');
    const countBadge = document.getElementById('candidate-count-badge');

    if (!tableBody) return;

    if (!candidates || candidates.length === 0) {
      tableBody.innerHTML = `
        <tr>
          <td colspan="8" class="text-center py-8 text-muted">
            No candidates ranked yet. Train a QSAR model first, then click "Rank Candidates".
          </td>
        </tr>
      `;
      if (countBadge) countBadge.textContent = '0 prioritized';
      return;
    }

    let filtered = candidates.filter(item => {
      const q = candidateTableState.searchQuery.toLowerCase();
      if (!q) return true;
      return (
        (item.compound_id && item.compound_id.toLowerCase().includes(q)) ||
        (item.compound_name && item.compound_name.toLowerCase().includes(q))
      );
    });

    const col = candidateTableState.sortColumn;
    const dir = candidateTableState.sortDirection === 'asc' ? 1 : -1;
    filtered.sort((a, b) => {
      let vA = a[col];
      let vB = b[col];
      return ((Number(vA) || 0) - (Number(vB) || 0)) * dir;
    });

    if (countBadge) {
      countBadge.textContent = `${filtered.length} candidates prioritized`;
    }

    const totalPages = Math.ceil(filtered.length / candidateTableState.pageSize) || 1;
    if (candidateTableState.page > totalPages) candidateTableState.page = totalPages;
    if (candidateTableState.page < 1) candidateTableState.page = 1;

    const startIdx = (candidateTableState.page - 1) * candidateTableState.pageSize;
    const pageItems = filtered.slice(startIdx, startIdx + candidateTableState.pageSize);

    tableBody.innerHTML = pageItems
      .map((item, idx) => {
        const rank = startIdx + idx + 1;
        const estIC50 = QSAR.pic50ToNanomolar(item.predicted_pic50);
        let badgeColor = '#ef4444';
        let badgeText = 'Predicted Inactive';
        if (item.predicted_pic50 >= 7.0) {
          badgeColor = '#10b981';
          badgeText = 'Predicted Active';
        } else if (item.predicted_pic50 >= 6.0) {
          badgeColor = '#f59e0b';
          badgeText = 'Predicted Intermediate';
        }

        const knownPic50 = item.actual_pic50 !== undefined ? item.actual_pic50.toFixed(2) : 'N/A';
        const delta =
          item.actual_pic50 !== undefined
            ? (item.predicted_pic50 - item.actual_pic50).toFixed(2)
            : '—';

        return `
        <tr>
          <td><strong class="rank-number">#${rank}</strong></td>
          <td><strong class="text-cyan">${Validator.sanitize(item.compound_id)}</strong></td>
          <td>${Validator.sanitize(item.compound_name)}</td>
          <td><strong class="text-teal">${item.predicted_pic50.toFixed(2)}</strong></td>
          <td>${estIC50 < 1000 ? estIC50.toFixed(1) + ' nM' : (estIC50 / 1000).toFixed(2) + ' µM'}</td>
          <td><span class="badge" style="background-color: ${badgeColor}22; color: ${badgeColor}; border: 1px solid ${badgeColor}44;">${badgeText}</span></td>
          <td>${knownPic50}</td>
          <td>${delta}</td>
        </tr>
      `;
      })
      .join('');

    if (paginationInfo) {
      paginationInfo.textContent = `Page ${candidateTableState.page} of ${totalPages}`;
    }
    if (prevBtn) prevBtn.disabled = candidateTableState.page <= 1;
    if (nextBtn) nextBtn.disabled = candidateTableState.page >= totalPages;
  }

  /**
   * Show Detailed Modal for a Compound
   * @param {string} compoundId
   */
  function showCompoundModal(compoundId) {
    const dataset = DataManager.getCurrentDataset();
    const compound = dataset.find(c => c.compound_id === compoundId);
    if (!compound) return;

    const modal = document.getElementById('compound-modal');
    const modalContent = document.getElementById('compound-modal-details');
    const canvas = document.getElementById('compound-structure-canvas');
    if (!modal || !modalContent) return;

    // Evaluate Lipinski violations
    let lipinskiViolations = 0;
    const lipinskiBadges = [];
    if (compound.molecular_weight > 500) {
      lipinskiViolations++;
      lipinskiBadges.push('MW > 500 Da');
    }
    if (compound.logP > 5.0) {
      lipinskiViolations++;
      lipinskiBadges.push('LogP > 5.0');
    }
    if (compound.hbd > 5) {
      lipinskiViolations++;
      lipinskiBadges.push('HBD > 5');
    }
    if (compound.hba > 10) {
      lipinskiViolations++;
      lipinskiBadges.push('HBA > 10');
    }

    const chemblLink = compound.compound_id.startsWith('CHEMBL')
      ? `https://www.ebi.ac.uk/chembl/compound_report_card/${compound.compound_id}/`
      : null;

    modalContent.innerHTML = `
      <div class="modal-header-info">
        <h3>${Validator.sanitize(compound.compound_name)} <span class="text-cyan">(${Validator.sanitize(compound.compound_id)})</span></h3>
        ${chemblLink ? `<a href="${chemblLink}" target="_blank" rel="noopener noreferrer" class="external-link">View in ChEMBL ↗</a>` : ''}
      </div>

      <div class="structure-render-box mb-4">
        <canvas id="modal-struct-canvas" width="480" height="240"></canvas>
      </div>

      <div class="info-grid">
        <div class="info-card">
          <span class="label">SMILES</span>
          <span class="value font-mono text-xs break-all">${Validator.sanitize(compound.smiles)}</span>
        </div>
        <div class="info-card">
          <span class="label">Experimental IC₅₀</span>
          <span class="value">${compound.activity_value} ${compound.activity_unit || 'nM'}</span>
        </div>
        <div class="info-card">
          <span class="label">Experimental pIC₅₀</span>
          <span class="value text-cyan font-bold">${compound.pic50.toFixed(2)}</span>
        </div>
        <div class="info-card">
          <span class="label">Activity Status</span>
          <span class="value">${compound.activity_class}</span>
        </div>
      </div>

      <h4 class="mt-4 mb-2 text-sm text-cyan uppercase tracking-wider">Molecular Descriptors</h4>
      <div class="descriptors-summary-table">
        <div class="desc-item"><span>Molecular Weight:</span> <strong>${compound.molecular_weight.toFixed(2)} Da</strong></div>
        <div class="desc-item"><span>Calculated LogP:</span> <strong>${compound.logP.toFixed(2)}</strong></div>
        <div class="desc-item"><span>H-Bond Donors (HBD):</span> <strong>${compound.hbd}</strong></div>
        <div class="desc-item"><span>H-Bond Acceptors (HBA):</span> <strong>${compound.hba}</strong></div>
        <div class="desc-item"><span>Rotatable Bonds:</span> <strong>${compound.rotatable_bonds}</strong></div>
        <div class="desc-item"><span>TPSA:</span> <strong>${compound.tpsa.toFixed(2)} Å²</strong></div>
        <div class="desc-item"><span>Total Rings:</span> <strong>${compound.ring_count}</strong></div>
        <div class="desc-item"><span>Aromatic Rings:</span> <strong>${compound.aromatic_ring_count}</strong></div>
        <div class="desc-item"><span>Fraction Csp3:</span> <strong>${compound.fraction_sp3.toFixed(2)}</strong></div>
      </div>

      <div class="lipinski-box mt-4 ${lipinskiViolations <= 1 ? 'lipinski-pass' : 'lipinski-fail'}">
        <strong>Lipinski Rule of 5:</strong> 
        ${lipinskiViolations === 0 ? '✓ 0 Violations (Drug-like oral bioavailability)' : `⚠️ ${lipinskiViolations} Violation(s) [${lipinskiBadges.join(', ')}]`}
      </div>

      <div class="disclaimer-note mt-3">
        <small class="text-muted">Bioactivity records retrieved from ChEMBL for target EGFR (ErbB1). Computational parameters require experimental biochemical validation.</small>
      </div>
    `;

    modal.classList.add('modal-open');

    // Draw 2D chemical structure in canvas
    setTimeout(() => {
      const structCanvas = document.getElementById('modal-struct-canvas');
      if (structCanvas) {
        Descriptors.renderMolecule(structCanvas, compound.smiles, { bgColor: '#0b1120' });
      }
    }, 50);
  }

  function closeModal() {
    const modal = document.getElementById('compound-modal');
    if (modal) modal.classList.remove('modal-open');
  }

  return {
    showToast,
    navigateTo,
    renderDatasetTable,
    renderCandidateTable,
    showCompoundModal,
    closeModal,
    getDatasetTableState: () => datasetTableState,
    getCandidateTableState: () => candidateTableState,
    getActiveSection: () => activeSectionId
  };
})();

// Attach to window if running in browser
if (typeof window !== 'undefined') {
  window.UI = UI;
}
