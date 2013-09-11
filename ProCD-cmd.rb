## Command-Line version of ProCD
## Written by Hanjo Kim
## 2013. 8. 9

## CHANGES
## 1. [2013. 9. 11] Added '-m suc' option. This option will mutate all lysines to negative charges.

require 'choice' # for option parsing
require File.dirname(__FILE__) + '/lib/procd'

PROGRAM_VERSION = "1.0.1"

Choice.options do
	banner "### ProCD command-line version #{PROGRAM_VERSION} ############################################\n"
	banner "###                                                                         ###\n"
	banner "###  Written by Hanjo Kim (hjkim@equisnzaroo.com)                           ###\n"
	banner "###  Only for use in EQUISnZAROO Inc.                                       ###\n"
	banner "###                                                                         ###\n"
	banner "###  Usage: ProCD-cmd [-iortmhv]                                            ###\n"
	banner "###  Please use --help option to get help message.                          ###\n"
	banner "###############################################################################"
	
	header ''
	header 'Specific options:'
	
	option :pdb_file, :required => true do
		short '-i'
		long '--input=PDB_FILE'
		desc 'Filename of protein in pdb format'
	end
	
	option :mode do
		short '-o'
		long '--mode=MODE'
		desc 'Program execution mode: NORMAL/AUTO_MUTATE'
		default "normal"
	end
	
	option :rsa_file do
		short '-r'
		long '--rsa=RSA_FILE'
		desc 'Filename of rsa file produced by naccess'
	end
	
	option :threshold do
		short '-t'
		long '--threshold==PERCENTAGE'
		desc 'Threshold of RSA'
		cast Float
		default 10.0
	end
	
	option :mutations do
		short '-m'
		long '--mutations=MUTATIONS'
		desc 'Specifies mutation sites, e.g. "D65A A109Q"'
    desc 'Multiple mutations should be double quoted.'
    desc '"suc" means "all Lys to negative"'
		validate /(suc)|([A-Z]{1}\d+[A-Z]{1}\s?)/
	end
	
	separator ''
	separator 'Common options: '
	
	option :help do
		short '-h'
		long '--help'
		desc 'Show this message'
	end
	
	option :version do
		short '-v'
		long '--version'
		desc 'Show version info'
		action do
			puts "ProCD-CMD.exe v#{PROGRAM_VERSION}"
		end
	end
end

pdb_file = Choice.choices[:pdb_file]
mode = Choice.choices[:mode]
rsa_file = Choice.choices[:rsa_file]
threshold = Choice.choices[:threshold]
mutations = Choice.choices[:mutations]

if mode.upcase == "NORMAL"
	output = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
	:threshold=>threshold, :mutations=> mutations)
	
	## Print the output to console
	puts "Result of ProCD calculation"
	puts "  protein file: #{pdb_file}"
	if rsa_file
		puts "  rsa file: #{rsa_file}"
		puts "  threshold: #{threshold.to_s}%"
	end
	if mutations
		puts "  mutation: #{mutations}"
	end
	puts "  executed at #{Time.now}"
	puts ""
	case output.fetch(:error_code)
		when 1
			puts "Protein file error!"
		when 2
			puts "Wrong mutation site!"
		when 0 # normal case
			puts "### Summary of Calculation ###"
			puts "  - Positive residues: #{output.fetch(:pos_size)}"
			puts "  - Negative Residues: #{output.fetch(:neg_size)}"
			puts "  - Magnitude of P sum vector: #{output.fetch(:pp_length)}"
			puts "  - Magnitude of N sum vector: #{output.fetch(:nn_length)}"
			puts "  - Angle between two vectors: #{output.fetch(:angle).round(2).to_s}"
			puts "  - ProCD Score: #{output.fetch(:score).round(3).to_s}"
	end
elsif mode.upcase == "AUTO_MUTATE"
	puts "Data preparing..."
	auto_mutate(:pdb_file => pdb_file, :rsa_file => rsa_file, :threshold=>threshold)
	puts "Result file written."
else
	puts "Please provide valid mode option."
	puts "To see the usage, 'ProCD-cmd --help'"
end