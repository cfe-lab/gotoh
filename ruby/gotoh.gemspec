Gem::Specification.new do |s|
  s.name = "gotoh"
  s.version = ENV['GOTOH_VERSION'] || '0.1.0.pre'
  s.summary = "CfE implementation of the Gotoh sequence alignment algorithm"
  s.files = [
    'lib/gotoh.rb',
    'ext/gotoh/gotoh.cpp'
  ]
  s.extensions = [
    'ext/gotoh/extconf.rb',
  ]
  s.authors = [
    "Conan Woods",
    "Jamie Kai",
    "David Rickett",
    "Richard Liang"
  ]
end
