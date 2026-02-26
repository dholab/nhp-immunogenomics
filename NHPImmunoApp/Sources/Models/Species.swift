import Foundation

/// Species metadata from alleles.json.
struct Species: Codable, Hashable, Sendable {
    let scientificName: String
    let commonName: String
    let taxon: Int
}
