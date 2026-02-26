// swift-tools-version: 6.2

import PackageDescription

let package = Package(
    name: "NHPImmunoApp",
    platforms: [.macOS(.v26)],
    targets: [
        .executableTarget(
            name: "NHPImmunoApp",
            path: "Sources"
        ),
        .testTarget(
            name: "NHPImmunoAppTests",
            dependencies: ["NHPImmunoApp"],
            path: "Tests"
        ),
    ]
)
