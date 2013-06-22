###############################################################################
# pI calc ver 0.02
# 2013. 06. 19.
# Written by Hanjo Kim

# Simple routine for calculation of pI value from protein sequence.
# Inspired by http://isoelectric.ovh.org/files/practise-isoelectric-point.html
# Requires BioRuby
#
# 1. Read the pdb file
# 2. Extract sequence using seqres method of BioRuby
# 3. Calculate the charge at certain pH
# 4. Change the pH until nq < E value (threshold)

# USAGE:
# calc_isoelectric_point pdb_file, chain_id

# Updates
# 1. exclude ssbonded cysteine residues from calculation

# TODO
# 1. Extract sequence from pdb atom list
# 2. Assign calculated pKa of each residue (from other file)
# 3. Exclude SSbonded cysteine residues from calculation **SOLVED**
###############################################################################

require 'bio'
require 'bio/db/pdb'

AminoAcids = Bio::AminoAcid.aa.keys[0..19].map {|r| Bio::AminoAcid.to_3(r)}

def get_seq_from_pdb pdb_file, chain_id
	protein = Bio::PDB.new(IO.read(pdb_file))
	protein.seqres(chain_id)
end

# res_name in 3-letter code
def count_res seq, res_name
	codes = seq.codes # codes.class = "Array"
	codes.count(res_name)
end

def percentage_res seq
  total = seq.length
  a = AminoAcids
  v = Hash.new
  a.each do |res|
    v[res] = ((count_res(seq, res) / total.to_f) * 100).round(2)
  end
  return v
end

def count_cys pdb_file, chain_id
	protein = Bio::PDB.new(IO.read(pdb_file))
	codes = protein.seqres(chain_id)
	total_cys_number = codes.count("Cys")
	all_ssbonded_cys = protein.record["SSBOND"]
	all_ssbonded_cys.keep_if {|c| (c.original_data[0][15] == chain_id)&&(c.original_data[0][29]==chain_id)}
	free_cys = total_cys_number - (all_ssbonded_cys.size * 2)
	return free_cys
end

def calc_nq ph, asp_number, glu_number, cys_number, tyr_number, his_number, lys_number, arg_number
	# calculate each charge
	qn1 = -1 / (1 + 10 ** (3.65 - ph)) 					# C-terminal charge
	qn2 = -1 * asp_number / (1 + 10 ** (3.9 -ph)) 		# Asp charge
	qn3 = -1 * glu_number / (1 + 10 ** (4.07 -ph)) 		# Glu charge
	qn4 = -1 * cys_number / (1 + 10 ** (8.18 - ph)) 	# Cys charge
	qn5 = -1 * tyr_number / (1 + 10 ** (10.46 - ph)) 	# Tyr charge
	
	qp1 = his_number / (1 + 10 ** (ph - 6.04)) 			# His charge
	qp2 = 1 / (1 + 10 ** (ph - 8.2)) 					# N-terminal charge
	qp3 = lys_number / (1 + 10 ** (ph - 10.54)) 		# Lys charge
	qp4 = arg_number / (1 + 10 ** (ph - 12.48)) 		# Arg charge
	
	# return nq value
	qn1 + qn2 + qn3 + qn4 + qn5 + qp1 + qp2 + qp3 + qp4
end

def calc_isoelectric_point pdb_file, chain_id
	# count residue numbers
	seq = get_seq_from_pdb pdb_file, chain_id
	seq_length = seq.length
	asp_number = count_res(seq, "Asp")
	glu_number = count_res(seq, "Glu")
	cys_number = count_cys(pdb_file, chain_id)
	tyr_number = count_res(seq, "Tyr")
	his_number = count_res(seq, "His")
	lys_number = count_res(seq, "Lys")
	arg_number = count_res(seq, "Arg")
	
	# initialize pH values
	ph = 6.5
	ph_prev = 0.0
	ph_next = 14.0
	e = 0.001 # nq value threshold. nq < e then quit the loop.
	
	nq = calc_nq(ph, asp_number, glu_number, cys_number, tyr_number, his_number, lys_number, arg_number)
	#puts "pH: 6.5 #{ph_prev.round(2)}~#{ph_next.round(2)} >> nq: #{nq.to_s}"
	
	until (nq.abs < e)
	#until (ph - ph_prev < e) && (ph_next - ph < e)
		if nq < 0
			temp = ph
			ph = ph - ((ph - ph_prev) / 2.0)
			ph_next = temp
			nq = calc_nq(ph, asp_number, glu_number, cys_number, tyr_number, his_number, lys_number, arg_number)
		else
			temp = ph
			ph = ph + ((ph_next - ph) / 2.0)
			ph_prev = temp
			nq = calc_nq(ph, asp_number, glu_number, cys_number, tyr_number, his_number, lys_number, arg_number)
		end
		#puts "pH: #{ph.round(4)} #{ph_prev.round(4)}~#{ph_next.round(4)} >> nq: #{nq.round(4).to_s}"
	end
	return ph.round(2)
end
