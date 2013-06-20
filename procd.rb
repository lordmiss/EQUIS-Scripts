require 'bio'
require 'bio/db/pdb'

def get_res_list_by_charge protein, opt
	list = Array.new
	neutrals = ["ALA", "GLY", "ILE", "LEU", "MET", "PRO", "VAL",
				"ASN", "CYS", "GLN", "SER", "THR", "PHE", "TRP",
				"TYR", "ARG", "HIS"]
	positives = ["ARG", "LYS"]
	negatives = ["ASP", "GLU"]
	protein.residues.each do |res|
		case opt.upcase
			when "P" then list << res if positives.include?(res.resName)
			when "N" then list << res if negatives.include?(res.resName)
			else list << res if neutrals.include?(res.resName)
		end
	end
	return list # returns an array of residues
end

def get_normalized_residue_vector protein, res
	p_vector = protein.centreOfGravity
	r_vector = res.centreOfGravity
	inner_vector = r_vector - p_vector
	inner_vector.normalize
end

def get_vector_sum protein, opt
	res_list = case opt.upcase
		when "P" then get_res_list_by_charge(protein, "p")
		when "N" then get_res_list_by_charge(protein, "n")
	end
	vectors = res_list.map {|res| get_normalized_residue_vector(protein, res)}
	vectors.inject {|sum, x| sum + x} # returns the sum of all vectors in the list
end

def filter_by_sa res_list, rsa_file, threshold
	list = res_list.map{|res| "RES #{res.resName} #{res.chain.id}#{sprintf('%4d', res.resSeq)}"}
	rsa = Array.new
	IO.foreach(rsa_file) do |line|
		list.each do |string|
			rsa << line[13..21].strip.to_f if line.start_with?(string)
		end
	end
	rsa.map! {|f| f < threshold ? true : false}
	(0..(rsa.size - 1)).each do |i|
		res_list.delete_at(i) if rsa[i]
	end
	return res_list
end