require 'bio'
require 'bio/db/pdb'

def get_protein pdb_file, chain_id
	protein = Bio::PDB.new(IO.read(pdb_file))
	protein.models[0][chain_id]
end

def get_atoms protein
	self.atoms
end

def get_distance a, b
	Math.sqrt((a.xyz.x - b.xyz.x)**2 + (a.xyz.y - b.xyz.y)**2 + (a.xyz.z - b.xyz.z)**2).round(3)
end

def get_cx_matrix atoms
	sphere = 4 / 3 * Math::PI * 1000
	size = atoms.size
	scores = Array.new(size, 0)
	matrix = Matrix.build(size, size) {|r, c| get_distance(atoms[r], atoms[c]) < 10.0 ? 1 : 0}.to_a
	(0..(size-1)).each do |i|
		value = ((sphere - ((matrix[i].count(1) -1) * 20.1)) / ((matrix[i].count(1) -1) * 20.1)).round(3)
		if value < 0
			scores[i] = 0.0
		else
			scores[i] = value
		end
	end
	return scores
end

def get_cx protein, a # 'a' should be "Bio::PDB::Record::ATOM" class
	atoms = get_atoms(protein)
	sphere = 4 / 3 * Math::PI * 1000
	int_atoms = 0
	atoms.each do |atom|
		distance = get_distance(a, atom)
		unless atom == a
			int_atoms = int_atoms + 1 if distance < 10.0
		end
	end
	return ((sphere - (int_atoms * 20.1)) / (int_atoms * 20.1)).round(3)
end

def get_cx_res protein, res # residue name: L111, C25, etc.
	atoms = get_atoms protein
	res_name = Bio::AminoAcid::Data::NAMES[res[0]].upcase
	res_seq = res[1..-1].to_i
	residue = protein.find_residue {|x| x.resName == res_name && x.resSeq == res_seq}[0]
	if residue == nil
		puts "Invalid residue name is provided. Please find the right name."
		return nil
	end
	res_atoms = protein.find_atom {|y| y.resName == res_name && y.resSeq == res_seq}
	sum = 0
	res_atoms.each do |atom|
		sum = sum + get_cx(protein, atom)
	end
	return (sum / res_atoms.size).round(3)
end

def get_cx_avg protein, cx_matrix
	cx_avg = Array.new(protein.residues.size, 0.0)
	resNo = 0
	atmNo = 0
	protein.residues.each do |residue|
		size = residue.atoms.size
		index = atmNo..(atmNo + size - 1)
		sum = 0.0
		index.each do |i|
			sum = sum + cx_matrix[i]
		end
		cx_avg[resNo] = (sum / size).round(3)
		resNo = resNo + 1
		atmNo = atmNo + size
	end
	return cx_avg
end