require 'simplecov'
require 'simplecov-cobertura'


SimpleCov.start do
  formatter SimpleCov::Formatter::MultiFormatter.new([
    SimpleCov::Formatter::HTMLFormatter, # Add HTML report for viewing/review.
    SimpleCov::Formatter::CoberturaFormatter # For CI
    ])
end

require 'minitest/autorun'
require 'minitest/reporters'

# if ENV['CI_COMMIT_REF']
#   run_id = ENV['CI_COMMIT_REF']
# elsif ENV["CI"]
if ENV['CI_COMMIT_REF']
  run_id = "#{ENV['CI_COMMIT_REF']} #{ENV['CI_PIPELINE_ID']}-#{ENV['CI_COMMIT_SHA']}"
else
  run_id = "LOCALBUILD"
end

Minitest::Reporters.use! [
  Minitest::Reporters::SpecReporter.new,
  Minitest::Reporters::JUnitReporter.new,
  Minitest::Reporters::HtmlReporter.new(
    :title => "cfe_gotoh Test Report #{run_id}",
    :erb_template => File.join(File.dirname(__FILE__), "templates/index.html.erb")
  )
]
