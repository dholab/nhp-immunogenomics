import Foundation

/// Manage provisional allele submissions via the GitHub API.
/// Creates branches, commits files, creates PRs, and triggers workflows.
actor ProvisionalService {
    private let api: GitHubAPI

    init(api: GitHubAPI) {
        self.api = api
    }

    /// Sanitize a string for use in a Git branch name.
    private func sanitizeBranchComponent(_ input: String) -> String {
        input.map { char -> String in
            switch char {
            case "*", ":", " ", "~", "^", "?", "[", "\\", ".":
                return "-"
            default:
                return String(char)
            }
        }.joined()
    }

    // MARK: - Read

    /// Fetch the current provisional manifest from main.
    func fetchManifest() async throws -> [ProvisionalAllele] {
        do {
            let (content, _) = try await api.getFileContent(path: "provisional/manifest.tsv")
            return TSVParser.parseManifest(content)
        } catch GitHubAPIError.httpError(404, _) {
            return []
        }
    }

    // MARK: - Add

    /// Submit new provisional alleles.
    /// Creates a branch, stages FASTA + config, creates a draft PR, triggers the workflow.
    func submitNewAlleles(species: String, locus: String, alleleClass: String,
                          seqType: String, submitter: String, notes: String,
                          fastaContent: String) async throws -> URL {
        let timestamp = Int(Date().timeIntervalSince1970)
        let branch = "provisional/\(sanitizeBranchComponent(species))-\(sanitizeBranchComponent(locus))-\(timestamp)"

        // Create branch from main
        let mainSHA = try await api.getBranchSHA("main")
        try await api.createBranch(name: branch, fromSHA: mainSHA)

        // Stage the config file
        let config: [String: String] = [
            "species": species, "locus": locus, "allele_class": alleleClass,
            "seq_type": seqType, "submitter": submitter, "notes": notes,
        ]
        let configData = try JSONEncoder().encode(config)
        try await api.commitFile(path: ".staging/config.json", content: configData,
                                  message: "Stage provisional submission config", branch: branch)

        // Stage the FASTA file
        let fastaData = Data(fastaContent.utf8)
        try await api.commitFile(path: ".staging/input.fasta", content: fastaData,
                                  message: "Stage provisional FASTA", branch: branch)

        // Create draft PR
        let records = FASTAParser.parse(fastaContent)
        let prBody = """
        ## Provisional Allele Submission

        | Field | Value |
        |-------|-------|
        | **Species** | \(species) |
        | **Locus** | \(locus) |
        | **Class** | \(alleleClass) |
        | **Seq type** | \(seqType) |
        | **Submitter** | \(submitter) |
        | **Sequences** | \(records.count) |
        | **Notes** | \(notes) |

        Submitted via NHP Immunogenomics app.
        Waiting for workflow to process and assign allele names.
        """

        let (_, prURL) = try await api.createPullRequest(
            title: "Add provisional: \(species)-\(locus) (\(records.count) sequence(s))",
            body: prBody,
            head: branch,
            draft: true
        )

        // Trigger the processing workflow on this branch
        try await api.triggerWorkflow(
            workflowFileName: "process-provisional-changes.yml",
            ref: branch
        )

        return prURL
    }

    // MARK: - Edit

    /// Edit a provisional allele's metadata (notes, submitter).
    func editAllele(_ updated: ProvisionalAllele) async throws -> URL {
        let timestamp = Int(Date().timeIntervalSince1970)
        let safeName = sanitizeBranchComponent(updated.name)
        let branch = "provisional/edit-\(safeName)-\(timestamp)"

        // Get current manifest
        let (content, manifestSHA) = try await api.getFileContent(path: "provisional/manifest.tsv")

        // Update the entry
        let newContent = TSVParser.updateInManifest(content, updated: updated)

        // Create branch and commit
        let mainSHA = try await api.getBranchSHA("main")
        try await api.createBranch(name: branch, fromSHA: mainSHA)
        try await api.commitFile(path: "provisional/manifest.tsv",
                                  content: Data(newContent.utf8),
                                  message: "Edit provisional allele \(updated.name)",
                                  branch: branch,
                                  existingSHA: manifestSHA)

        // Create PR
        let (_, prURL) = try await api.createPullRequest(
            title: "Edit provisional: \(updated.name)",
            body: "Updated metadata for provisional allele \(updated.name).\n\nSubmitted via NHP Immunogenomics app.",
            head: branch
        )

        // Trigger rebuild
        try await api.triggerWorkflow(
            workflowFileName: "process-provisional-changes.yml",
            ref: branch
        )

        return prURL
    }

    // MARK: - Delete

    /// Delete provisional alleles by name.
    func deleteAlleles(_ names: Set<String>) async throws -> URL {
        let timestamp = Int(Date().timeIntervalSince1970)
        let label = names.count == 1 ? names.first! : "\(names.count)-alleles"
        let safeName = sanitizeBranchComponent(label)
        let branch = "provisional/delete-\(safeName)-\(timestamp)"

        // Get current manifest
        let (content, manifestSHA) = try await api.getFileContent(path: "provisional/manifest.tsv")
        let alleles = TSVParser.parseManifest(content)

        // Identify affected FASTA files
        let toDelete = alleles.filter { names.contains($0.name) }
        let affectedFiles = Set(toDelete.map(\.sequenceFile))

        // Create branch
        let mainSHA = try await api.getBranchSHA("main")
        try await api.createBranch(name: branch, fromSHA: mainSHA)

        // Commit updated manifest
        let newManifest = TSVParser.removeFromManifest(content, names: names)
        try await api.commitFile(path: "provisional/manifest.tsv",
                                  content: Data(newManifest.utf8),
                                  message: "Remove \(names.count) provisional allele(s)",
                                  branch: branch,
                                  existingSHA: manifestSHA)

        // Update affected FASTA files (remove deleted sequences)
        for file in affectedFiles {
            let fastaPath = "provisional/sequences/\(file)"
            do {
                let (fastaContent, fastaSHA) = try await api.getFileContent(path: fastaPath, ref: branch)
                let records = FASTAParser.parse(fastaContent)
                let namesToRemove = Set(toDelete.filter { $0.sequenceFile == file }.map(\.name))
                let kept = records.filter { !namesToRemove.contains($0.name) }
                if kept.isEmpty {
                    try await api.deleteFile(path: fastaPath,
                                              message: "Remove empty FASTA \(file)",
                                              branch: branch, sha: fastaSHA)
                } else {
                    let newFasta = FASTAParser.format(kept)
                    try await api.commitFile(path: fastaPath,
                                              content: Data(newFasta.utf8),
                                              message: "Remove sequences from \(file)",
                                              branch: branch, existingSHA: fastaSHA)
                }
            } catch GitHubAPIError.httpError(404, _) {
                continue
            }
        }

        // Create PR
        let nameList = names.sorted().joined(separator: ", ")
        let (_, prURL) = try await api.createPullRequest(
            title: "Delete provisional: \(label)",
            body: "Removing provisional allele(s): \(nameList)\n\nSubmitted via NHP Immunogenomics app.",
            head: branch
        )

        // Trigger rebuild
        try await api.triggerWorkflow(
            workflowFileName: "process-provisional-changes.yml",
            ref: branch
        )

        return prURL
    }
}
