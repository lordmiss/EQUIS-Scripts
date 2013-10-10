require 'bio'
$:.unshift File.join(File.dirname(__FILE__), 'lib')
require 'calc_pi'
require 'procd'
require 'dpx'

def get_protein_name_and_chain_id name
  values = Hash.new
  base = File.basename(name, ".pdb").split("_")
  values["name"] = base[0]
  values["chain_id"] = base[1]
  return values
end

pdbs = Dir["../../Temp/cytokine_pdb_chain/*.pdb"]
rsas = Dir["../../Temp/cytokine_pdb_chain/*.rsa"]
header = ["file_name", "pdb_code", "chain_id", "pi", "ProCD_no_filter", "ProCD_10", "ProCD_dpx"]
output = Array.new

pdbs.each do |pdb|
  protein = move_pdb_to_origin(pdb)
  make_25_proteins(protein)
  make_dpx_pdb_file(pdb)
end

pdbs.each do |pdb|
  rsa = rsas[pdbs.index(pdb)]
  name = get_protein_name_and_chain_id(pdb)["name"]
  chain_id = get_protein_name_and_chain_id(pdb)["chain_id"]
  line = Array.new
  line << File.basename(pdb)
  line << name
  line << chain_id
  line << calc_isoelectric_point(pdb, chain_id)
  line << calc_procd_score(:pdb_file=>pdb)[:score]
  line << calc_procd_score(:pdb_file=>pdb, :rsa_file=>rsa)[:score]
  line << calc_procd_score(:pdb_file=>pdb, :dpx=>true)
  output << line
  print "."
end

CSV.open("total.csv", "wb") do |csv|
  csv << header
  output.each {|l| csv << l}
end