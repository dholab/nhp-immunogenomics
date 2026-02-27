import Testing
@testable import NHPImmunoApp

@Test func parseSingleRecord() {
    let fasta = ">seq1\nATGCATGC\nGCTAGCTA\n"
    let records = FASTAParser.parse(fasta)
    #expect(records.count == 1)
    #expect(records[0].name == "seq1")
    #expect(records[0].sequence == "ATGCATGCGCTAGCTA")
}

@Test func parseMultipleRecords() {
    let fasta = """
    >allele1
    ATGC
    >allele2
    GCTA
    """
    let records = FASTAParser.parse(fasta)
    #expect(records.count == 2)
    #expect(records[0].name == "allele1")
    #expect(records[1].name == "allele2")
}

@Test func parseCRLFLineEndings() {
    let fasta = ">seq1\r\nATGC\r\nGCTA\r\n"
    let records = FASTAParser.parse(fasta)
    #expect(records.count == 1)
    #expect(records[0].name == "seq1")
    #expect(records[0].sequence == "ATGCGCTA")
}

@Test func validateGoodSequence() {
    #expect(FASTAParser.validateSequence("ATGCNRYSWKMBDHV") == nil)
}

@Test func validateBadSequence() {
    let error = FASTAParser.validateSequence("ATGX")
    #expect(error != nil)
    #expect(error!.contains("X"))
}

@Test func validateEmptySequence() {
    let error = FASTAParser.validateSequence("")
    #expect(error != nil)
}

@Test func formatRoundtrip() {
    let records = [FASTARecord(name: "test", sequence: String(repeating: "ATGC", count: 20))]
    let formatted = FASTAParser.format(records)
    let parsed = FASTAParser.parse(formatted)
    #expect(parsed.count == 1)
    #expect(parsed[0].name == "test")
    #expect(parsed[0].sequence == records[0].sequence)
}

@Test func formatTrailingNewline() {
    let records = [FASTARecord(name: "seq", sequence: "ATGC")]
    let formatted = FASTAParser.format(records)
    #expect(formatted.hasSuffix("\n"))
}
