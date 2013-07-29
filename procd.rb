require 'bio'
require 'bio/db/pdb'
require 'csv'

# Edit!

# TODO
# 1. Chain ID parsing

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
	opts = {:threshold => 10.0, :rsa_file=>nil, :mutation=>nil }.merge(opts) # default threshold is 10.0     <- default option. if no user input, use default parameter
	pdb_file = opts.fetch(:pdb_file) # name of pdb file  user
	rsa_file = opts.fetch(:rsa_file) # name of rsa file  optional
	threshold = opts.fetch(:threshold) # RSA threshold	  optional
	mutation = opts.fetch(:mutation) # mutation	  optional
	
	# parse protein
	begin
		protein = Bio::PDB.new(IO.read(pdb_file))
	rescue
		puts "Protein file error!"
		return nil
	end
	residues = protein.residues
	
	if rsa_file == nil
		pos = get_res_list_by_charge(:protein=>protein, :charge=>"P")
		neg = get_res_list_by_charge(:protein=>protein, :charge=>"N")
	else
		pos = get_res_list_by_charge(:protein=>protein, :rsa_file=>rsa_file, :threshold=>threshold, :charge=>"P")
		neg = get_res_list_by_charge(:protein=>protein, :rsa_file=>rsa_file, :threshold=>threshold, :charge=>"N")
	end
	
	pos_size = pos.size
	neg_size = neg.size
	
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
			pos_size = pos_size - 1
			neg_size = neg_size + 1
		elsif POSITIVES.include?(orig) && NEUTRALS.include?(mut) # positive delete
			pp = pp - m_vec
			pos_size = pos_size - 1
		elsif NEGATIVES.include?(orig) && POSITIVES.include?(mut) # negative -> positive
			pp = pp + m_vec
			nn = nn - m_vec
			pos_size = pos_size + 1
			neg_size = neg_size - 1
		elsif NEGATIVES.include?(orig) && NEUTRALS.include?(mut) # negative delete
			nn = nn - m_vec
			neg_size = neg_size - 1
		elsif NEUTRALS.include?(orig) && POSITIVES.include?(mut) # positive add
			pp = pp + m_vec
			pos_size = pos_size + 1
		elsif NEUTRALS.include?(orig) && NEGATIVES.include?(mut) # negative add
			nn = nn + m_vec
			neg_size = neg_size + 1
		end
	end
	
	# Result output
	puts "Positive Residues: #{pos_size}"
	puts "Negative Residues: #{neg_size}"
	puts "Magnitude of P sum vector: #{pp.magnitude.round(3).to_s}"
	puts "Magnitude of N sum vector: #{nn.magnitude.round(3).to_s}"
	puts "Angle between two vectors: #{angle.round(2).to_s}"
	puts "ProCD Score: #{(pp.magnitude/nn.magnitude).round(3)}"
	
	values = {:pos_size => pos_size, :neg_size => neg_size,
			:pp_length => pp.magnitude.round(3),
			:nn_length => nn.magnitude.round(3),
			:angle => angle.round(2),
			:score => (pp.magnitude/nn.magnitude).round(3)}
	return values
end

def calc_angle a, b
	Math::acos((a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (a.magnitude * b.magnitude)) * 180 / Math::PI
end

def auto_mutate(opts={})
	# Option parsing
	opts = {:threshold => 10.0, :rsa_file=>nil}.merge(opts)
	pdb_file = opts.fetch(:pdb_file) # name of pdb file 
	chain_id = opts.fetch(:chain_id) # chain ID
	rsa_file = opts.fetch(:rsa_file) # name of rsa file
	threshold = opts.fetch(:threshold).to_f # RSA threshold

	rsa_option = case threshold
		when nil
			"filter"
		else
			"none"
	end
	
	time_stamp = Time.new.strftime("%Y%m%d_%H%M%S_%L")
	output_file_name = "#{pdb_file}_#{time_stamp}.txt"
	begin
		protein = Bio::PDB.new(IO.read(pdb_file)).models[0][chain_id]
	rescue
		puts "Protein file error!"
		return nil
	end
	residues = protein.residues
	
	# calculate original procd score
	orig = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file, :threshold => threshold)
	
	# prepare output text
	output = [[time_stamp],
			["RSA option : #{rsa_option}"],
			["Naccess : #{threshold.to_s}"],
			["Naccess class : AllAtomREL"],
			["default", "default", "#{orig.fetch(:pos_size)}", "#{orig.fetch(:neg_size)}",
			 "#{orig.fetch(:pp_length)}", "#{orig.fetch(:nn_length)}", "#{orig.fetch(:score)}"]
			]
	# mutate the protein twice per residue
	residues.each do |res|
		id = res.id
		res_name3 = res.resName.capitalize
		res_name = Bio::AminoAcid.to_1(res_name3) # 1-letter amino acid code
				
		pos_mutation = "#{res_name}#{id.to_s}K"
		neg_mutation = "#{res_name}#{id.to_s}D"
		del_mutation = "#{res_name}#{id.to_s}A"
		if NEUTRALS.include?(res.resName)
			pos = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file,
				:threshold => threshold, :mutation=>pos_mutation)
			neg = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file,
				:threshold => threshold, :mutation=>neg_mutation)
			output << ["[#{res_name3}]#{id.to_s}:#{chain_id}", "positive", "#{pos.fetch(:pos_size)}",
			"#{pos.fetch(:neg_size)}", "#{pos.fetch(:pp_length)}", "#{pos.fetch(:nn_length)}", "#{pos.fetch(:score)}"]
			output << ["[#{res_name3}]#{id.to_s}:#{chain_id}", "negative", "#{neg.fetch(:pos_size)}",
			"#{neg.fetch(:neg_size)}", "#{neg.fetch(:pp_length)}", "#{neg.fetch(:nn_length)}", "#{neg.fetch(:score)}"]
		elsif POSITIVES.include?(res.resName)
			del = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file,
				:threshold => threshold, :mutation=>del_mutation)
			neg = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file,
				:threshold => threshold, :mutation=>neg_mutation)
			output << ["[#{res_name3}]#{id.to_s}:#{chain_id}", "delete", "#{del.fetch(:pos_size)}",
			"#{del.fetch(:neg_size)}", "#{del.fetch(:pp_length)}", "#{del.fetch(:nn_length)}", "#{del.fetch(:score)}"]
			output << ["[#{res_name3}]#{id.to_s}:#{chain_id}", "pos->neg", "#{neg.fetch(:pos_size)}",
			"#{neg.fetch(:neg_size)}", "#{neg.fetch(:pp_length)}", "#{neg.fetch(:nn_length)}", "#{neg.fetch(:score)}"]
		else
			del = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file,
				:threshold => threshold, :mutation=>del_mutation)
			pos = calc_procd_score(:pdb_file => pdb_file, rsa_file => rsa_file,
				:threshold => threshold, :mutation=>pos_mutation)
			output << ["[#{res_name3}]#{id.to_s}:#{chain_id}", "delete", "#{del.fetch(:pos_size)}",
			"#{del.fetch(:neg_size)}", "#{del.fetch(:pp_length)}", "#{del.fetch(:nn_length)}", "#{del.fetch(:score)}"]
			output << ["[#{res_name3}]#{id.to_s}:#{chain_id}", "neg->pos", "#{pos.fetch(:pos_size)}",
			"#{pos.fetch(:neg_size)}", "#{pos.fetch(:pp_length)}", "#{pos.fetch(:nn_length)}", "#{pos.fetch(:score)}"]
		end
	end
	
	CSV.open(output_file_name, "wb") do |csv|
		output.each do |line|
			csv << line
		end
	end
end