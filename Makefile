.PHONY: all build build-release build-production check test doc clean prepublish

# Default target
all: build test doc

build:
	@cargo build
	@cargo clippy

check:
	@cargo check
	@cargo clippy

clean:
	@cargo clean

build-release:
	@cargo build --release

build-production:
	@cargo build --profile production

test:
	@cargo test

doc:
	@cargo doc --no-deps

prepublish:
	@cargo publish --dry-run
	@cargo package --list
	@echo "see https://doc.rust-lang.org/cargo/reference/publishing.html"
