import Foundation

/// Parse and write tab-separated values files (manifest.tsv format).
enum TSVParser {

    /// Parse TSV content into ProvisionalAllele records.
    /// Skips the header row and any blank lines.
    static func parseManifest(_ content: String) -> [ProvisionalAllele] {
        let normalized = content.normalizedLineEndings
        let lines = normalized.split(separator: "\n", omittingEmptySubsequences: true)
            .map(String.init)
        guard lines.count > 1 else { return [] }
        return lines.dropFirst().compactMap { ProvisionalAllele(tsvRow: $0) }
    }

    /// Serialize provisional alleles back to TSV with header row.
    static func formatManifest(_ alleles: [ProvisionalAllele]) -> String {
        let header = ProvisionalAllele.columns.joined(separator: "\t")
        let rows = alleles.map(\.tsvRow)
        return ([header] + rows).joined(separator: "\n") + "\n"
    }

    /// Remove specific alleles from a manifest string by name.
    static func removeFromManifest(_ content: String, names: Set<String>) -> String {
        let alleles = parseManifest(content).filter { !names.contains($0.name) }
        return formatManifest(alleles)
    }

    /// Update a specific allele in a manifest string.
    static func updateInManifest(_ content: String, updated: ProvisionalAllele) -> String {
        var alleles = parseManifest(content)
        if let idx = alleles.firstIndex(where: { $0.name == updated.name }) {
            alleles[idx] = updated
        }
        return formatManifest(alleles)
    }
}
