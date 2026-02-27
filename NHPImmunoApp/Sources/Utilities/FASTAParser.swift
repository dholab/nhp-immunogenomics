import Foundation

/// Normalize line endings from CRLF or CR to LF.
extension String {
    var normalizedLineEndings: String {
        replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
    }
}

/// A single FASTA record: header name and nucleotide sequence.
struct FASTARecord: Sendable, Hashable {
    let name: String
    let sequence: String
}

/// Lightweight FASTA parser for local validation before upload.
enum FASTAParser {

    /// Valid IUPAC nucleotide characters (case-insensitive).
    private static let validNucleotides: Set<Character> = Set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")

    /// Parse FASTA content into records.
    static func parse(_ content: String) -> [FASTARecord] {
        let normalized = content.normalizedLineEndings
        var records: [FASTARecord] = []
        var currentName: String?
        var currentSeq: [String] = []

        for line in normalized.split(separator: "\n", omittingEmptySubsequences: false) {
            let trimmed = line.trimmingCharacters(in: .whitespaces)
            if trimmed.hasPrefix(">") {
                if let name = currentName {
                    records.append(FASTARecord(name: name, sequence: currentSeq.joined()))
                }
                currentName = String(trimmed.dropFirst()).trimmingCharacters(in: .whitespaces)
                currentSeq = []
            } else if !trimmed.isEmpty {
                currentSeq.append(trimmed)
            }
        }
        if let name = currentName {
            records.append(FASTARecord(name: name, sequence: currentSeq.joined()))
        }
        return records
    }

    /// Validate a sequence string. Returns nil if valid, or an error message.
    static func validateSequence(_ sequence: String) -> String? {
        if sequence.isEmpty {
            return "Sequence is empty"
        }
        let invalid = sequence.filter { !validNucleotides.contains($0) }
        if !invalid.isEmpty {
            let chars = Set(invalid).sorted().prefix(5).map { String($0) }.joined(separator: ", ")
            return "Invalid nucleotide characters: \(chars)"
        }
        return nil
    }

    /// Validate all records. Returns a list of error messages (empty = valid).
    static func validateRecords(_ records: [FASTARecord]) -> [String] {
        var errors: [String] = []
        if records.isEmpty {
            errors.append("No sequences found in FASTA input")
            return errors
        }
        for record in records {
            if record.name.isEmpty {
                errors.append("Record has empty header")
            }
            if let err = validateSequence(record.sequence) {
                errors.append("\(record.name): \(err)")
            }
        }
        return errors
    }

    /// Format records back into FASTA text with 70-character line wrapping.
    /// Uses incremental index advancement to avoid O(n²) string indexing.
    static func format(_ records: [FASTARecord]) -> String {
        var output = ""
        for record in records {
            output += ">\(record.name)\n"
            var index = record.sequence.startIndex
            let end = record.sequence.endIndex
            while index < end {
                let lineEnd = record.sequence.index(index, offsetBy: 70, limitedBy: end) ?? end
                output += record.sequence[index..<lineEnd]
                output += "\n"
                index = lineEnd
            }
        }
        return output
    }
}
