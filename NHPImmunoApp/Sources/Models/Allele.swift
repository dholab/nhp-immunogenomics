import Foundation

/// A single allele record from alleles.json.
/// Uses CodingKeys to map the minified JSON keys to readable Swift property names.
struct Allele: Identifiable, Codable, Hashable, Sendable {
    let accession: String
    let name: String
    let locus: String
    let alleleClass: String
    let species: String
    let project: String
    let dateAssigned: String
    let dateModified: String
    var previousNames: [String]?
    var provisional: Bool?
    var submitter: String?
    var seqType: String?

    var id: String { accession }

    enum CodingKeys: String, CodingKey {
        case accession = "a"
        case name = "n"
        case locus = "l"
        case alleleClass = "c"
        case species = "s"
        case project = "p"
        case dateAssigned = "da"
        case dateModified = "dm"
        case previousNames = "prev"
        case provisional = "prov"
        case submitter = "sub"
        case seqType = "st"
    }

    var projectLabel: String {
        switch project {
        case "MHC": "MHC"
        case "NHKIR": "KIR"
        default: project.capitalized
        }
    }

    var isProvisional: Bool { provisional == true }
}
