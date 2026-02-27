import Foundation

/// A provisional allele parsed from manifest.tsv.
/// Column order: name, species, locus, class, seq_type, sequence_file, submitter, date_added, notes
struct ProvisionalAllele: Identifiable, Hashable, Sendable {
    var name: String
    var species: String
    var locus: String
    var alleleClass: String
    var seqType: String
    var sequenceFile: String
    var submitter: String
    var dateAdded: String
    var notes: String

    var id: String { name }

    /// TSV column names in manifest order.
    static let columns = ["name", "species", "locus", "class", "seq_type",
                           "sequence_file", "submitter", "date_added", "notes"]

    /// Serialize to a TSV row (tab-separated values).
    var tsvRow: String {
        [name, species, locus, alleleClass, seqType,
         sequenceFile, submitter, dateAdded, notes].joined(separator: "\t")
    }
}

extension ProvisionalAllele {
    /// Parse from a TSV row (tab-separated values).
    /// Handles notes fields that may contain tab characters by joining remaining fields.
    init?(tsvRow: String) {
        // Only trim newlines, not whitespace — tabs are meaningful TSV delimiters
        let trimmed = tsvRow.trimmingCharacters(in: .newlines)
        guard !trimmed.isEmpty else { return nil }
        let fields = trimmed.split(separator: "\t", omittingEmptySubsequences: false)
            .map(String.init)
        guard fields.count >= 9 else { return nil }
        self.name = fields[0]
        self.species = fields[1]
        self.locus = fields[2]
        self.alleleClass = fields[3]
        self.seqType = fields[4]
        self.sequenceFile = fields[5]
        self.submitter = fields[6]
        self.dateAdded = fields[7]
        self.notes = fields[8...].joined(separator: "\t")
    }
}
