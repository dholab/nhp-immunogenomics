import Foundation

/// Top-level structure of alleles.json served from GitHub Pages.
struct AlleleIndex: Codable, Sendable {
    let generated: String
    let mhcVersion: String
    let nhkirVersion: String
    let species: [String: Species]
    let alleles: [Allele]
    let loci: [String: [String]]

    enum CodingKeys: String, CodingKey {
        case generated
        case mhcVersion = "mhc_version"
        case nhkirVersion = "nhkir_version"
        case species, alleles, loci
    }
}
