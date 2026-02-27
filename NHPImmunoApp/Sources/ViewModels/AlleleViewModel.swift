import Foundation
import Observation

@MainActor
@Observable
final class AlleleViewModel {
    // Data
    var alleles: [Allele] = []
    var speciesMap: [String: Species] = [:]
    var lociByProject: [String: [String]] = [:]
    var mhcVersion = ""
    var nhkirVersion = ""
    var generated = ""

    // Filters
    var searchText = ""
    var selectedSpecies: String?
    var selectedLocus: String?
    var selectedClass: String?
    var selectedProject: String?
    var showProvisionalOnly = false

    // State
    var isLoading = false
    var errorMessage: String?
    var selectedAlleleID: Allele.ID?

    // Cached results — updated via applyFilters()
    private(set) var filteredAlleles: [Allele] = []
    private(set) var availableSpecies: [String] = []
    private(set) var availableLoci: [String] = []

    /// Persistent service instance to preserve cache across loads.
    private var alleleService: AlleleService?
    private var serviceOwner: String?
    private var serviceRepo: String?

    /// O(1) lookup for selected allele.
    private var alleleLookup: [Allele.ID: Allele] = [:]

    /// The currently selected allele, resolved from the ID.
    var selectedAllele: Allele? {
        guard let id = selectedAlleleID else { return nil }
        return alleleLookup[id]
    }

    var alleleCount: String {
        let filtered = filteredAlleles.count
        let total = alleles.count
        return filtered == total ? "\(total) alleles" : "\(filtered) of \(total) alleles"
    }

    /// Recompute filteredAlleles, availableSpecies, and availableLoci from current filter state.
    /// Call this after changing any filter property or after loading new data.
    func applyFilters() {
        var result = alleles
        if !searchText.isEmpty {
            let query = searchText.lowercased()
            result = result.filter { allele in
                allele.name.lowercased().contains(query) ||
                allele.accession.lowercased().contains(query) ||
                (allele.previousNames?.contains(where: { $0.lowercased().contains(query) }) ?? false)
            }
        }
        if let sp = selectedSpecies { result = result.filter { $0.species == sp } }
        if let lo = selectedLocus { result = result.filter { $0.locus == lo } }
        if let cl = selectedClass { result = result.filter { $0.alleleClass == cl } }
        if let pr = selectedProject { result = result.filter { $0.project == pr } }
        if showProvisionalOnly { result = result.filter(\.isProvisional) }
        filteredAlleles = result

        availableSpecies = Set(alleles.map(\.species)).sorted()

        var locusSource = alleles
        if let sp = selectedSpecies { locusSource = locusSource.filter { $0.species == sp } }
        if let proj = selectedProject { locusSource = locusSource.filter { $0.project == proj } }
        availableLoci = Set(locusSource.map(\.locus)).sorted()
    }

    func load(owner: String, repo: String) async {
        isLoading = true
        errorMessage = nil
        defer { isLoading = false }

        if alleleService == nil || serviceOwner != owner || serviceRepo != repo {
            alleleService = AlleleService(owner: owner, repo: repo)
            serviceOwner = owner
            serviceRepo = repo
        }
        guard let service = alleleService else { return }
        do {
            let index = try await service.fetchIndex(forceRefresh: true)
            alleles = index.alleles
            alleleLookup = Dictionary(alleles.map { ($0.id, $0) }, uniquingKeysWith: { _, b in b })
            speciesMap = index.species
            lociByProject = index.loci
            mhcVersion = index.mhcVersion
            nhkirVersion = index.nhkirVersion
            generated = index.generated
            applyFilters()
            // Release the service cache — data is now in the ViewModel
            await service.clearCache()
        } catch {
            errorMessage = error.localizedDescription
        }
    }

    var hasActiveFilters: Bool {
        selectedSpecies != nil || selectedLocus != nil
            || selectedClass != nil || selectedProject != nil
            || showProvisionalOnly
    }

    func clearFilters() {
        searchText = ""
        selectedSpecies = nil
        selectedLocus = nil
        selectedClass = nil
        selectedProject = nil
        showProvisionalOnly = false
        applyFilters()
    }
}
