/* =========================================================
   NHP Immunogenomics Database - Application Logic
   Vanilla JS, no dependencies, no build step.
   ========================================================= */

(function () {
  "use strict";

  // -------------------------------------------------------
  // Constants
  // -------------------------------------------------------
  var ROWS_PER_PAGE = 100;
  var DATA_URL = "alleles.json";

  // -------------------------------------------------------
  // State
  // -------------------------------------------------------
  var state = {
    alleles: [],
    species: {},
    loci: {},
    metadata: {},
    filtered: [],
    sorted: [],
    selected: new Set(),
    sortKey: "n",
    sortAsc: true,
    page: 0,
    filters: {
      database: "All",
      species: "",
      mhcClass: "",
      locus: "",
      search: "",
    },
  };

  // -------------------------------------------------------
  // DOM references (cached after DOMContentLoaded)
  // -------------------------------------------------------
  var dom = {};

  // -------------------------------------------------------
  // Initialization
  // -------------------------------------------------------
  document.addEventListener("DOMContentLoaded", function () {
    cacheDom();
    bindEvents();
    loadData();
  });

  function cacheDom() {
    dom.loadingState = document.getElementById("loading-state");
    dom.errorState = document.getElementById("error-state");
    dom.errorDetail = document.getElementById("error-detail");
    dom.retryBtn = document.getElementById("retry-btn");
    dom.appContent = document.getElementById("app-content");
    dom.versionInfo = document.getElementById("version-info");

    dom.databaseToggle = document.querySelectorAll(".toggle-btn[data-value]");
    dom.speciesSelect = document.getElementById("species-select");
    dom.classSelect = document.getElementById("class-select");
    dom.locusSelect = document.getElementById("locus-select");
    dom.nameSearch = document.getElementById("name-search");
    dom.clearFiltersBtn = document.getElementById("clear-filters-btn");

    dom.resultsCount = document.getElementById("results-count");
    dom.downloadSection = document.getElementById("download-section");
    dom.downloadLabel = document.getElementById("download-label");
    dom.downloadBtns = document.querySelectorAll(".btn-download");
    dom.selectAll = document.getElementById("select-all");
    dom.resultsBody = document.getElementById("results-body");
    dom.emptyState = document.getElementById("empty-state");
    dom.pagination = document.getElementById("pagination");
    dom.pageInfo = document.getElementById("page-info");
    dom.prevPage = document.getElementById("prev-page");
    dom.nextPage = document.getElementById("next-page");
    dom.sortBtns = document.querySelectorAll(".sort-btn");
  }

  function bindEvents() {
    // Database toggle
    dom.databaseToggle.forEach(function (btn) {
      btn.addEventListener("click", function () {
        setDatabase(btn.dataset.value);
      });
    });

    // Select filters
    dom.speciesSelect.addEventListener("change", function () {
      state.filters.species = this.value;
      onFilterChange();
    });

    dom.classSelect.addEventListener("change", function () {
      state.filters.mhcClass = this.value;
      onFilterChange();
    });

    dom.locusSelect.addEventListener("change", function () {
      state.filters.locus = this.value;
      onFilterChange();
    });

    // Search with debounce
    var searchTimer = null;
    dom.nameSearch.addEventListener("input", function () {
      var val = this.value;
      clearTimeout(searchTimer);
      searchTimer = setTimeout(function () {
        state.filters.search = val;
        onFilterChange();
      }, 200);
    });

    // Clear filters
    dom.clearFiltersBtn.addEventListener("click", clearAllFilters);

    // Retry
    dom.retryBtn.addEventListener("click", loadData);

    // Usage "More..." toggle
    var usageLink = document.getElementById("usage-more-link");
    var usageDetail = document.getElementById("usage-detail");
    if (usageLink && usageDetail) {
      usageLink.addEventListener("click", function (e) {
        e.preventDefault();
        var showing = !usageDetail.hidden;
        usageDetail.hidden = showing;
        usageLink.textContent = showing ? "More\u2026" : "Less";
      });
    }

    // Sort
    dom.sortBtns.forEach(function (btn) {
      btn.addEventListener("click", function () {
        var key = btn.dataset.sort;
        if (state.sortKey === key) {
          state.sortAsc = !state.sortAsc;
        } else {
          state.sortKey = key;
          state.sortAsc = true;
        }
        applySorting();
        renderTable();
        updateSortIndicators();
      });
    });

    // Pagination
    dom.prevPage.addEventListener("click", function () {
      if (state.page > 0) {
        state.page--;
        renderTable();
        updatePagination();
        scrollToTable();
      }
    });

    dom.nextPage.addEventListener("click", function () {
      var maxPage = Math.ceil(state.filtered.length / ROWS_PER_PAGE) - 1;
      if (state.page < maxPage) {
        state.page++;
        renderTable();
        updatePagination();
        scrollToTable();
      }
    });

    // Select all checkbox
    dom.selectAll.addEventListener("change", function () {
      var start = state.page * ROWS_PER_PAGE;
      var end = Math.min(start + ROWS_PER_PAGE, state.filtered.length);
      var pageData = state.filtered.slice(start, end);

      if (dom.selectAll.checked) {
        pageData.forEach(function (a) { state.selected.add(a.a); });
      } else {
        pageData.forEach(function (a) { state.selected.delete(a.a); });
      }
      renderTable();
      updateSelectionUI();
    });

    // Download buttons
    dom.downloadBtns.forEach(function (btn) {
      btn.addEventListener("click", function () {
        handleDownload(btn.dataset.format, btn);
      });
    });
  }

  // -------------------------------------------------------
  // Data Loading
  // -------------------------------------------------------
  function loadData() {
    dom.loadingState.hidden = false;
    dom.errorState.hidden = true;
    dom.appContent.hidden = true;

    fetch(DATA_URL)
      .then(function (res) {
        if (!res.ok) {
          throw new Error("HTTP " + res.status + ": " + res.statusText);
        }
        return res.json();
      })
      .then(function (data) {
        state.alleles = data.alleles || [];
        state.species = data.species || {};
        state.loci = data.loci || {};
        state.metadata = {
          generated: data.generated || "",
          mhcVersion: data.mhc_version || "",
          nhkirVersion: data.nhkir_version || "",
        };
        initializeUI();
      })
      .catch(function (err) {
        dom.loadingState.hidden = true;
        dom.errorState.hidden = false;
        dom.errorDetail.textContent = err.message;
      });
  }

  function initializeUI() {
    dom.loadingState.hidden = true;
    dom.appContent.hidden = false;

    // Version info
    var generated = state.metadata.generated
      ? new Date(state.metadata.generated).toLocaleDateString("en-US", {
          year: "numeric",
          month: "long",
          day: "numeric",
        })
      : "unknown";
    dom.versionInfo.textContent =
      "IPD-MHC v" +
      state.metadata.mhcVersion +
      " | IPD-NHKIR v" +
      state.metadata.nhkirVersion +
      " | Generated " +
      generated;

    populateSpeciesDropdown();
    populateLocusDropdown();
    updateClassSelectState();
    onFilterChange();
    updateSortIndicators();
  }

  // -------------------------------------------------------
  // Dropdown Population
  // -------------------------------------------------------
  function populateSpeciesDropdown() {
    // Clear existing options beyond the default
    dom.speciesSelect.length = 1;

    var codes = Object.keys(state.species).sort();
    codes.forEach(function (code) {
      var sp = state.species[code];
      var opt = document.createElement("option");
      opt.value = code;
      opt.textContent = code + " - " + sp.commonName;
      dom.speciesSelect.appendChild(opt);
    });
  }

  function populateLocusDropdown() {
    var currentValue = state.filters.locus;
    dom.locusSelect.length = 1;

    var db = state.filters.database;
    var lociSet = new Set();

    if (db === "All" || db === "MHC") {
      (state.loci.MHC || []).forEach(function (l) {
        lociSet.add(l);
      });
    }
    if (db === "All" || db === "NHKIR") {
      (state.loci.NHKIR || []).forEach(function (l) {
        lociSet.add(l);
      });
    }

    var lociArr = Array.from(lociSet).sort();
    lociArr.forEach(function (locus) {
      var opt = document.createElement("option");
      opt.value = locus;
      opt.textContent = locus;
      dom.locusSelect.appendChild(opt);
    });

    // Restore selection if still valid
    if (lociArr.indexOf(currentValue) !== -1) {
      dom.locusSelect.value = currentValue;
    } else {
      dom.locusSelect.value = "";
      state.filters.locus = "";
    }
  }

  function updateClassSelectState() {
    var db = state.filters.database;
    // Disable MHC Class dropdown when only NHKIR is selected
    var disabled = db === "NHKIR";
    dom.classSelect.disabled = disabled;
    if (disabled) {
      dom.classSelect.value = "";
      state.filters.mhcClass = "";
    }
  }

  // -------------------------------------------------------
  // Filter Logic
  // -------------------------------------------------------
  function setDatabase(value) {
    state.filters.database = value;

    // Update toggle UI
    dom.databaseToggle.forEach(function (btn) {
      var isActive = btn.dataset.value === value;
      btn.classList.toggle("active", isActive);
      btn.setAttribute("aria-checked", isActive ? "true" : "false");
    });

    updateClassSelectState();
    populateLocusDropdown();
    onFilterChange();
  }

  function clearAllFilters() {
    state.filters.database = "All";
    state.filters.species = "";
    state.filters.mhcClass = "";
    state.filters.locus = "";
    state.filters.search = "";

    dom.databaseToggle.forEach(function (btn) {
      var isAll = btn.dataset.value === "All";
      btn.classList.toggle("active", isAll);
      btn.setAttribute("aria-checked", isAll ? "true" : "false");
    });

    dom.speciesSelect.value = "";
    dom.classSelect.value = "";
    dom.classSelect.disabled = false;
    dom.nameSearch.value = "";

    populateLocusDropdown();
    onFilterChange();
  }

  function onFilterChange() {
    applyFilters();
    applySorting();
    state.page = 0;
    state.selected.clear();
    renderTable();
    updatePagination();
    updateResultsCount();
    updateSelectionUI();
  }

  function applyFilters() {
    var db = state.filters.database;
    var sp = state.filters.species;
    var cls = state.filters.mhcClass;
    var loc = state.filters.locus;
    var search = state.filters.search.toLowerCase().trim();

    state.filtered = state.alleles.filter(function (allele) {
      // Database
      if (db !== "All" && allele.p !== db) return false;

      // Species
      if (sp && allele.s !== sp) return false;

      // MHC Class
      if (cls && allele.c !== cls) return false;

      // Locus
      if (loc && allele.l !== loc) return false;

      // Search: match name or previous designations
      if (search) {
        var nameMatch = allele.n.toLowerCase().indexOf(search) !== -1;
        if (!nameMatch) {
          var prevMatch = false;
          if (allele.prev) {
            for (var i = 0; i < allele.prev.length; i++) {
              if (allele.prev[i].toLowerCase().indexOf(search) !== -1) {
                prevMatch = true;
                break;
              }
            }
          }
          if (!prevMatch) return false;
        }
      }

      return true;
    });
  }

  function applySorting() {
    var key = state.sortKey;
    var asc = state.sortAsc;

    state.filtered.sort(function (a, b) {
      var va = a[key] || "";
      var vb = b[key] || "";

      if (typeof va === "string") {
        va = va.toLowerCase();
        vb = vb.toLowerCase();
      }

      if (va < vb) return asc ? -1 : 1;
      if (va > vb) return asc ? 1 : -1;
      return 0;
    });
  }

  // -------------------------------------------------------
  // Rendering
  // -------------------------------------------------------
  function renderTable() {
    var start = state.page * ROWS_PER_PAGE;
    var end = Math.min(start + ROWS_PER_PAGE, state.filtered.length);
    var pageData = state.filtered.slice(start, end);

    if (state.filtered.length === 0) {
      dom.resultsBody.innerHTML = "";
      dom.emptyState.hidden = false;
      dom.pagination.hidden = true;
      return;
    }

    dom.emptyState.hidden = true;
    dom.pagination.hidden = state.filtered.length <= ROWS_PER_PAGE;

    // Update select-all checkbox state
    updateSelectAllCheckbox(pageData);

    // Build rows using document fragment for performance
    var fragment = document.createDocumentFragment();

    for (var i = 0; i < pageData.length; i++) {
      var allele = pageData[i];
      var tr = document.createElement("tr");
      if (state.selected.has(allele.a)) tr.classList.add("row-selected");

      // Checkbox
      var tdCheck = document.createElement("td");
      tdCheck.className = "cell-select";
      var cb = document.createElement("input");
      cb.type = "checkbox";
      cb.checked = state.selected.has(allele.a);
      cb.setAttribute("aria-label", "Select " + allele.n);
      cb.dataset.acc = allele.a;
      cb.addEventListener("change", (function (acc) {
        return function (e) {
          if (e.target.checked) {
            state.selected.add(acc);
          } else {
            state.selected.delete(acc);
          }
          e.target.closest("tr").classList.toggle("row-selected", e.target.checked);
          updateSelectAllCheckbox(pageData);
          updateSelectionUI();
        };
      })(allele.a));
      tdCheck.appendChild(cb);
      tr.appendChild(tdCheck);

      // Accession
      var tdAcc = document.createElement("td");
      tdAcc.className = "cell-accession";
      tdAcc.textContent = allele.a;
      tr.appendChild(tdAcc);

      // Name + previous designations
      var tdName = document.createElement("td");
      var nameSpan = document.createElement("span");
      nameSpan.className = "cell-name";
      nameSpan.textContent = allele.n;
      tdName.appendChild(nameSpan);

      if (allele.prev && allele.prev.length > 0) {
        var prevSpan = document.createElement("span");
        prevSpan.className = "cell-prev";
        prevSpan.textContent = "prev: " + allele.prev.join(", ");
        tdName.appendChild(prevSpan);
      }
      tr.appendChild(tdName);

      // Locus
      var tdLocus = document.createElement("td");
      tdLocus.textContent = allele.l;
      tr.appendChild(tdLocus);

      // Class
      var tdClass = document.createElement("td");
      tdClass.className = "cell-class";
      tdClass.textContent = allele.c;
      tr.appendChild(tdClass);

      // Species
      var tdSpecies = document.createElement("td");
      tdSpecies.textContent = allele.s;
      if (state.species[allele.s]) {
        tdSpecies.title = state.species[allele.s].scientificName;
      }
      tr.appendChild(tdSpecies);

      // Project
      var tdProject = document.createElement("td");
      tdProject.textContent = allele.p;
      tr.appendChild(tdProject);

      // Date Modified
      var tdDate = document.createElement("td");
      tdDate.className = "cell-date";
      tdDate.textContent = allele.dm || "";
      tr.appendChild(tdDate);

      fragment.appendChild(tr);
    }

    dom.resultsBody.innerHTML = "";
    dom.resultsBody.appendChild(fragment);
  }

  function updatePagination() {
    var total = state.filtered.length;
    var start = state.page * ROWS_PER_PAGE + 1;
    var end = Math.min((state.page + 1) * ROWS_PER_PAGE, total);
    var maxPage = Math.max(0, Math.ceil(total / ROWS_PER_PAGE) - 1);

    if (total === 0) {
      dom.pagination.hidden = true;
      return;
    }

    dom.pagination.hidden = total <= ROWS_PER_PAGE;
    dom.pageInfo.textContent =
      "Page " + (state.page + 1) + " of " + (maxPage + 1);
    dom.prevPage.disabled = state.page === 0;
    dom.nextPage.disabled = state.page >= maxPage;
  }

  function updateResultsCount() {
    var total = state.filtered.length;
    if (total === 0) {
      dom.resultsCount.textContent = "No results";
    } else {
      var start = state.page * ROWS_PER_PAGE + 1;
      var end = Math.min((state.page + 1) * ROWS_PER_PAGE, total);
      dom.resultsCount.textContent =
        "Showing " +
        start.toLocaleString() +
        "\u2013" +
        end.toLocaleString() +
        " of " +
        total.toLocaleString() +
        " results";
    }
  }

  function updateSortIndicators() {
    dom.sortBtns.forEach(function (btn) {
      if (btn.dataset.sort === state.sortKey) {
        btn.setAttribute(
          "aria-sort",
          state.sortAsc ? "ascending" : "descending"
        );
      } else {
        btn.removeAttribute("aria-sort");
      }
    });
  }

  function scrollToTable() {
    var el = document.querySelector(".results-table-wrapper");
    if (el) {
      el.scrollIntoView({ behavior: "smooth", block: "start" });
    }
  }

  // -------------------------------------------------------
  // Selection helpers
  // -------------------------------------------------------
  function updateSelectAllCheckbox(pageData) {
    if (!pageData || pageData.length === 0) {
      dom.selectAll.checked = false;
      dom.selectAll.indeterminate = false;
      return;
    }
    var checkedCount = 0;
    pageData.forEach(function (a) {
      if (state.selected.has(a.a)) checkedCount++;
    });
    dom.selectAll.checked = checkedCount === pageData.length;
    dom.selectAll.indeterminate = checkedCount > 0 && checkedCount < pageData.length;
  }

  function updateSelectionUI() {
    var count = state.selected.size;
    if (count > 0) {
      dom.downloadLabel.innerHTML =
        "Download <strong>" + count + " selected</strong>:" +
        ' <button type="button" class="btn-clear-selection" id="clear-selection-btn">clear</button>';
      document.getElementById("clear-selection-btn").addEventListener("click", function () {
        state.selected.clear();
        renderTable();
        updateSelectionUI();
      });
    } else {
      dom.downloadLabel.textContent = "Download filtered:";
    }
  }

  /** Return the alleles to use for download: selected if any, otherwise all filtered. */
  function getDownloadAlleles() {
    if (state.selected.size > 0) {
      return state.alleles.filter(function (a) { return state.selected.has(a.a); });
    }
    return state.filtered;
  }

  // -------------------------------------------------------
  // Downloads - all derived from GenBank files
  // -------------------------------------------------------
  function handleDownload(format, btn) {
    if (btn.classList.contains("is-loading")) return;
    var dlAlleles = getDownloadAlleles();
    if (dlAlleles.length === 0) return;

    if (format === "csv") {
      downloadCSV();
      return;
    }

    setButtonLoading(btn, true);

    // Group download alleles by (project, species, locus)
    var groups = groupAlleles(dlAlleles);

    // Fetch .gb files for each group
    var fetches = groups.map(function (g) {
      var url = "data/" + g.project + "/" + g.species + "/" + g.locus + ".gb";
      return fetchFile(url).then(function (text) {
        return { text: text, accessions: g.accessions };
      });
    });

    Promise.all(fetches)
      .then(function (results) {
        var output = "";

        results.forEach(function (r) {
          if (!r.text) return;
          var gbRecords = splitGenBankRecords(r.text, r.accessions);
          if (format === "genbank") {
            output += gbRecords.map(function (rec) { return rec.raw; }).join("");
          } else if (format === "fasta") {
            output += gbRecordsToFasta(gbRecords);
          } else if (format === "protein") {
            output += gbRecordsToProteinFasta(gbRecords);
          } else if (format === "coding") {
            output += gbRecordsToCodingFasta(gbRecords);
          }
        });

        if (!output.trim()) {
          alert("No sequence data found for the current selection in this format.");
          return;
        }

        var filename = buildDownloadFilename(format);

        downloadBlob(output, filename, "text/plain");
      })
      .catch(function (err) {
        console.error("Download failed:", err);
        alert("Download failed: " + err.message);
      })
      .finally(function () {
        setButtonLoading(btn, false);
      });
  }

  function groupAlleles(alleles) {
    var map = {};
    alleles.forEach(function (a) {
      var key = a.p + "/" + a.s + "/" + a.l;
      if (!map[key]) {
        map[key] = { project: a.p, species: a.s, locus: a.l, accessions: new Set() };
      }
      map[key].accessions.add(a.a);
    });
    return Object.values(map);
  }

  function fetchFile(url) {
    return fetch(url)
      .then(function (res) {
        if (!res.ok) {
          if (res.status === 404) return "";
          throw new Error("HTTP " + res.status + " fetching " + url);
        }
        return res.text();
      })
      .catch(function (err) {
        console.warn("Could not fetch " + url + ":", err.message);
        return "";
      });
  }

  // -------------------------------------------------------
  // GenBank parsing helpers
  // -------------------------------------------------------

  /**
   * Split a multi-record GenBank file, filter by accession set.
   * Returns array of {accession, name, raw, sequence, translation, cdsRanges}.
   */
  function splitGenBankRecords(text, accessionSet) {
    if (!text) return [];
    var rawRecords = text.split(/\/\/\n?/);
    var results = [];

    rawRecords.forEach(function (block) {
      var trimmed = block.trim();
      if (!trimmed) return;

      var accMatch = trimmed.match(/^ACCESSION\s+(\S+)/m);
      if (!accMatch || !accessionSet.has(accMatch[1])) return;

      var acc = accMatch[1];
      var defMatch = trimmed.match(/^DEFINITION\s+(.+)/m);
      var defLine = defMatch ? defMatch[1].replace(/\.$/, "").trim() : acc;
      // Extract just the allele name (before the comma)
      var name = defLine.indexOf(",") !== -1 ? defLine.substring(0, defLine.indexOf(",")).trim() : defLine;

      results.push({
        accession: acc,
        name: name,
        raw: trimmed + "\n//\n",
        sequence: extractOriginSequence(trimmed),
        translation: extractTranslation(trimmed),
        cdsRanges: extractCdsRanges(trimmed),
      });
    });

    return results;
  }

  /** Extract the nucleotide sequence from the ORIGIN block. */
  function extractOriginSequence(record) {
    var idx = record.indexOf("\nORIGIN");
    if (idx === -1) return "";
    var block = record.substring(idx);
    // Remove "ORIGIN" header line, line numbers, and spaces; uppercase
    return block.replace(/^[^\n]*\n/, "").replace(/[\s0-9\/]/g, "").toUpperCase();
  }

  /** Extract /translation="..." from CDS feature. */
  function extractTranslation(record) {
    var match = record.match(/\/translation="([^"]+)"/);
    if (!match) return "";
    return match[1].replace(/\s+/g, "");
  }

  /** Extract CDS join coordinates as array of [start, end] (0-based). */
  function extractCdsRanges(record) {
    // Match CDS line and its join(...) or complement(join(...)) or simple span
    var cdsMatch = record.match(/^\s{5}CDS\s+(.*)/m);
    if (!cdsMatch) return null;

    // Collect continuation lines (indented further than qualifier lines)
    var lines = record.split("\n");
    var cdsLineIdx = -1;
    for (var i = 0; i < lines.length; i++) {
      if (/^\s{5}CDS\s+/.test(lines[i])) { cdsLineIdx = i; break; }
    }
    if (cdsLineIdx === -1) return null;

    var locStr = lines[cdsLineIdx].replace(/^\s{5}CDS\s+/, "");
    for (var j = cdsLineIdx + 1; j < lines.length; j++) {
      var line = lines[j];
      // Stop at next qualifier or feature
      if (/^\s{21}\//.test(line) || /^\s{5}\S/.test(line) || /^[A-Z]/.test(line)) break;
      locStr += line.trim();
    }

    // Parse ranges from join(1..73,204..473,...) or simple 1..1098
    var ranges = [];
    var rangePattern = /(\d+)\.\.(\d+)/g;
    var m;
    while ((m = rangePattern.exec(locStr)) !== null) {
      ranges.push([parseInt(m[1], 10) - 1, parseInt(m[2], 10)]); // Convert to 0-based start
    }
    return ranges.length > 0 ? ranges : null;
  }

  /** Build a FASTA header: >{name} IPD:{accession} */
  function fastaHeader(r) {
    return ">" + r.name + " IPD:" + r.accession;
  }

  /** Convert parsed GenBank records to nucleotide FASTA (full sequence). */
  function gbRecordsToFasta(records) {
    var out = "";
    records.forEach(function (r) {
      if (!r.sequence) return;
      out += fastaHeader(r) + "\n";
      out += wrapSequence(r.sequence) + "\n";
    });
    return out;
  }

  /** Convert parsed GenBank records to protein FASTA from CDS /translation. */
  function gbRecordsToProteinFasta(records) {
    var out = "";
    records.forEach(function (r) {
      if (!r.translation) return;
      out += fastaHeader(r) + "\n";
      out += wrapSequence(r.translation) + "\n";
    });
    return out;
  }

  /** Convert parsed GenBank records to coding FASTA by extracting CDS regions. */
  function gbRecordsToCodingFasta(records) {
    var out = "";
    records.forEach(function (r) {
      if (!r.sequence || !r.cdsRanges) return;
      var coding = "";
      r.cdsRanges.forEach(function (range) {
        coding += r.sequence.substring(range[0], range[1]);
      });
      if (!coding) return;
      out += fastaHeader(r) + "\n";
      out += wrapSequence(coding) + "\n";
    });
    return out;
  }

  /** Wrap a sequence string to 70 characters per line. */
  function wrapSequence(seq) {
    var lines = [];
    for (var i = 0; i < seq.length; i += 70) {
      lines.push(seq.substring(i, i + 70));
    }
    return lines.join("\n");
  }

  /** Build a descriptive filename based on active filters and format. */
  function buildDownloadFilename(format) {
    var parts = [];

    // Database
    var db = state.filters.database;
    if (db === "All") {
      parts.push("IPD-MHC_NHKIR");
    } else if (db === "MHC") {
      parts.push("IPD-MHC");
    } else {
      parts.push("IPD-NHKIR");
    }

    // Species
    if (state.filters.species) parts.push(state.filters.species);

    // MHC Class
    if (state.filters.mhcClass) parts.push("class" + state.filters.mhcClass);

    // Locus
    if (state.filters.locus) parts.push(state.filters.locus);

    // Search term (sanitized)
    if (state.filters.search) {
      var sanitized = state.filters.search.trim().replace(/[^a-zA-Z0-9*:-]/g, "_").substring(0, 30);
      if (sanitized) parts.push(sanitized);
    }

    var base = parts.join("_");

    // Format-specific suffix
    var formatNames = {
      genbank: ".gb",
      fasta: "_genomic.fasta",
      coding: "_cds.fasta",
      protein: "_protein.fasta",
      csv: "_metadata.csv",
    };

    return base + (formatNames[format] || ".txt");
  }

  function downloadCSV() {
    var headers = [
      "Accession",
      "Name",
      "Locus",
      "Class",
      "Species",
      "Scientific Name",
      "Common Name",
      "Project",
      "Date Assigned",
      "Date Modified",
      "Previous Designations",
    ];

    var rows = [headers.join(",")];
    var dlAlleles = getDownloadAlleles();

    dlAlleles.forEach(function (a) {
      var sp = state.species[a.s] || {};
      var prev = a.prev ? a.prev.join("; ") : "";
      var row = [
        csvEscape(a.a),
        csvEscape(a.n),
        csvEscape(a.l),
        csvEscape(a.c),
        csvEscape(a.s),
        csvEscape(sp.scientificName || ""),
        csvEscape(sp.commonName || ""),
        csvEscape(a.p),
        csvEscape(a.da || ""),
        csvEscape(a.dm || ""),
        csvEscape(prev),
      ];
      rows.push(row.join(","));
    });

    var filename = buildDownloadFilename("csv");

    downloadBlob(rows.join("\n") + "\n", filename, "text/csv");
  }

  function csvEscape(value) {
    if (typeof value !== "string") value = String(value);
    if (value.indexOf(",") !== -1 || value.indexOf('"') !== -1 || value.indexOf("\n") !== -1) {
      return '"' + value.replace(/"/g, '""') + '"';
    }
    return value;
  }

  function downloadBlob(text, filename, mimeType) {
    var blob = new Blob([text], { type: mimeType + ";charset=utf-8" });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.style.display = "none";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  function setButtonLoading(btn, loading) {
    if (loading) {
      btn.classList.add("is-loading");
      btn.setAttribute("aria-disabled", "true");
    } else {
      btn.classList.remove("is-loading");
      btn.removeAttribute("aria-disabled");
    }
  }
})();
