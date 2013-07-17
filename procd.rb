require 'bio'
require 'bio/db/pdb'

NEUTRALS = ["ALA", "GLY", "ILE", "LEU", "MET", "PRO", "VAL",
			"ASN", "CYS", "GLN", "SER", "THR", "PHE", "TRP",
			"TYR", "ARG", "HIS"]
POSITIVES = ["ARG", "LYS"]
NEGATIVES = ["ASP", "GLU"]

def get_res_list_by_charge(opts = {})
	# Option parsing
	opts = {:threshold => 10.0, :rsa_file => nil }.merge(opts)
	protein = opts.fetch(:protein)
	charge = opts.fetch(:charge)
	rsa_file = opts.fetch(:rsa_file)
	threshold = opts.fetch(:threshold)
	
	list = Array.new
	neutrals = NEUTRALS
	positives = POSITIVES
	negatives = NEGATIVES
	
	residues = protein.residues.sort_by {|x| x.resSeq}
	if (rsa_file != nil)
		rsa = parse_rsa_file rsa_file
		del = Array.new
		rsa.each do |l|
			del << l[0] if l[1] < threshold
		end
		residues.delete_if {|res| del.include?(res.resSeq)}
	end
	
	residues.each do |res|
		case charge.upcase
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

def get_vector_sum protein, res_list
	vectors = res_list.map {|res| get_normalized_residue_vector(protein, res)}
	vectors.inject {|sum, x| sum + x} # returns the sum of all vectors in the list
end

def parse_rsa_file rsa_file
	values = Array.new
	IO.foreach(rsa_file) do |line|
		items = line.split(" ")
		if items[0] == "RES"
			values << [items[3].to_i, items[5].to_f]
		end
	end
	a = values.sort_by {|x| x[0]}
	return a
end

def calc_procd_score(opts = {})
	# option parsing
	opts = {:threshold => 10.0, :rsa_file=>nil, :mutation=>nil }.merge(opts) # default threshold is 10.0
	pdb_file = opts.fetch(:pdb_file) # name of pdb file 
	rsa_file = opts.fetch(:rsa_file) # name of rsa file
	threshold = opts.fetch(:threshold) # RSA threshold
	mutation = opts.fetch(:mutation) # mutation
	
	# parse protein
	protein = Bio::PDB.new(IO.read(pdb_file))
	residues = protein.residues
	
	if rsa_file == nil
		pos = get_res_list_by_charge(:protein=>protein, :charge=>"P")
		neg = get_res_list_by_charge(:protein=>protein, :charge=>"N")
	else
		pos = get_res_list_by_charge(:protein=>protein, :rsa_file=>rsa_file, :threshold=>threshold, :charge=>"P")
		neg = get_res_list_by_charge(:protein=>protein, :rsa_file=>rsa_file, :threshold=>threshold, :charge=>"N")
	end
	
	pp = get_vector_sum(protein, pos)
	nn = get_vector_sum(protein, neg)
	angle = calc_angle pp, nn

	# mutation
	if mutation != nil
		# parse mutation, mutation should have this form, E421Q
		orig = Bio::AminoAcid::Data::NAMES[mutation[0]].upcase # 3-letter code
		mut = Bio::AminoAcid::Data::NAMES[mutation[-1]].upcase # 3-letter code
		site = mutation[1..-2].to_i # mutation site number
		
		res = protein.find_residue{|res| res.resName == orig && res.resSeq == site}[0]
		m_vec = get_normalized_residue_vector protein, res
		if POSITIVES.include?(orig) && NEGATIVES.include?(mut) # positive -> negative
			pp = pp - m_vec
			nn = nn + m_vec
		elsif POSITIVES.include?(orig) && NEUTRALS.include?(mut) # positive delete
			pp = pp - m_vec
		elsif NEGATIVES.include?(orig) && POSITIVES.include?(mut) # negative -> positive
			pp = pp + m_vec
			nn = nn - m_vec
		elsif NEGATIVES.include?(orig) && NEUTRALS.include?(mut) # negative delete
			nn = nn - m_vec
		elsif NEUTRALS.include?(orig) && POSITIVES.include?(mut) # positive add
			pp = pp + m_vec
		elsif NEUTRALS.include?(orig) && NEGATIVES.include?(mut) # negative add
			nn = nn + m_vec
		end
	end
	
	# Result output
	puts "Positive Residues: #{pos.size}"
	puts "Negative Residues: #{neg.size}"
	puts "Magnitude of P sum vector: #{pp.magnitude.round(3).to_s}"
	puts "Magnitude of N sum vector: #{nn.magnitude.round(3).to_s}"
	puts "Angle between two vectors: #{angle.round(2).to_s}"
	puts "ProCD Score: #{(pp.magnitude/nn.magnitude).round(3)}"
end

def calc_angle a, b
	Math::acos((a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (a.magnitude * b.magnitude)) * 180 / Math::PI
end