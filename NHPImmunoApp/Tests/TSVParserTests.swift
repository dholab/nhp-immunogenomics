import Testing
@testable import NHPImmunoApp

@Test func parseManifest() {
    let tsv = """
    name\tspecies\tlocus\tclass\tseq_type\tsequence_file\tsubmitter\tdate_added\tnotes
    Mamu-A1*026:new01\tMamu\tA1\tI\tcoding\tMamu/A1.fasta\tJ. Karl\t2026-01-15\ttest note
    """
    let alleles = TSVParser.parseManifest(tsv)
    #expect(alleles.count == 1)
    #expect(alleles[0].name == "Mamu-A1*026:new01")
    #expect(alleles[0].species == "Mamu")
    #expect(alleles[0].locus == "A1")
    #expect(alleles[0].submitter == "J. Karl")
    #expect(alleles[0].notes == "test note")
}

@Test func parseManifestCRLF() {
    let tsv = "name\tspecies\tlocus\tclass\tseq_type\tsequence_file\tsubmitter\tdate_added\tnotes\r\nA\tsp\tL\tI\tcoding\tf\ts\t2026-01-01\tnote\r\n"
    let alleles = TSVParser.parseManifest(tsv)
    #expect(alleles.count == 1)
    #expect(alleles[0].name == "A")
    #expect(alleles[0].notes == "note")
}

@Test func formatManifest() {
    let allele = ProvisionalAllele(
        name: "Test-A*001:new01", species: "Test", locus: "A",
        alleleClass: "I", seqType: "coding", sequenceFile: "Test/A.fasta",
        submitter: "Tester", dateAdded: "2026-01-01", notes: ""
    )
    let formatted = TSVParser.formatManifest([allele])
    let parsed = TSVParser.parseManifest(formatted)
    #expect(parsed.count == 1)
    #expect(parsed[0].name == allele.name)
}

@Test func removeFromManifest() {
    let tsv = """
    name\tspecies\tlocus\tclass\tseq_type\tsequence_file\tsubmitter\tdate_added\tnotes
    A\tsp\tL\tI\tcoding\tf\ts\t2026-01-01\t
    B\tsp\tL\tI\tcoding\tf\ts\t2026-01-01\t
    C\tsp\tL\tI\tcoding\tf\ts\t2026-01-01\t
    """
    let result = TSVParser.removeFromManifest(tsv, names: ["B"])
    let remaining = TSVParser.parseManifest(result)
    #expect(remaining.count == 2)
    #expect(remaining.map(\.name) == ["A", "C"])
}
