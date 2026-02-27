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

    /// Persistent service instance to preserve cache across loads.
    private var alleleService: AlleleService?

    /// The currently selected allele, resolved from the ID.
    var selectedAllele: Allele? {
        guard let id = selectedAlleleID else { return nil }
        return alleles.first { $0.id == id }
    }

    var availableSpecies: [String] {
        Set(alleles.map(\.species)).sorted()
    }

    var availableLoci: [String] {
        var filtered = alleles
        if let sp = selectedSpecies { filtered = filtered.filter { $0.species == sp } }
        if let proj = selectedProject { filtered = filtered.filter { $0.project == proj } }
        return Set(filtered.map(\.locus)).sorted()
    }

    var filteredAlleles: [Allele] {
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
        return result
    }

    var alleleCount: String {
        let filtered = filteredAlleles.count
        let total = alleles.count
        return filtered == total ? "\(total) alleles" : "\(filtered) of \(total) alleles"
    }

    func load(owner: String, repo: String) async {
        isLoading = true
        errorMessage = nil
        defer { isLoading = false }

        if alleleService == nil {
            alleleService = AlleleService(owner: owner, repo: repo)
        }
        guard let service = alleleService else { return }
        do {
            let index = try await service.fetchIndex(forceRefresh: true)
            alleles = index.alleles
            speciesMap = index.species
            lociByProject = index.loci
            mhcVersion = index.mhcVersion
            nhkirVersion = index.nhkirVersion
            generated = index.generated
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
    }
}
