Gem::Specification.new do |s|
  s.name = "cfe_gotoh"
  s.version = ENV['CFE_GOTOH_VERSION'] || '0.1.0.pre'
  s.summary = "CfE implementation of the Gotoh sequence alignment algorithm"
  s.files = [
    'lib/cfe_gotoh.rb',
    'ext/cfe_gotoh/cfe_gotoh.cpp'
  ]
  s.extensions = [
    'ext/cfe_gotoh/extconf.rb',
  ]
  s.authors = [
    "Conan Woods",
    "Jamie Kai",
    "David Rickett",
    "Richard Liang"
  ]
  # Associate this gem with a GitHub repo; this allows the gems to be
  # automatically associated with the repo when pushed to GitHub Packages.
  s.metadata = {
    "github_repo" => "ssh://github.com/cfe-lab/gotoh"
  }
end
