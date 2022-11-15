//
//  SHA2.swift
//
//  Created by Karl-Johan Alm on 2022-10-25.
//

import Foundation

final class SHA2Hash {
    private var inner: sha256 = sha256()
    var data: Data {
        get { withUnsafeBytes(of: inner) { Data($0) } }
        set { inner = newValue.withUnsafeBytes { $0.load(as: sha256.self) } }
    }
    init(data input: Data) {
        data = input
    }
    init(finalizeContext: inout sha256_ctx) {
        sha256_done(&inner, &finalizeContext)
    }
}

final class SHA2Context {
    private var initialized = true
    private var inner: sha256_ctx = sha256_ctx()
    init() {
        sha256_init(&inner)
    }
    func reinitialize() -> Void {
        assert(!initialized)
        initialized = true
        sha256_init(&inner)
    }
    func update(data: Data) -> Void {
        assert(initialized)
        sha256_update(&inner, [UInt8](data), data.count)
    }
    func done() -> Data {
        assert(initialized)
        initialized = false
        return SHA2Hash(finalizeContext: &inner).data
    }
    var ready: Bool { initialized }
}

// TODO: sha256_double64 wrapper, if needed
// TODO: sha256_midstate wrapper, if needed
