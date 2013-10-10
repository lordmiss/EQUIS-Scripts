#require 'bio'
require 'bio/db/pdb'
require 'csv'

# TODO
# 1. Chain ID parsing
# 2. fix the auto_mutate method (it does not use rsa_file option) : SOLVED

# CHANGES
# 1. [2013. 9. 11] Added '-m suc' option. This option will mutate all lysines to negative charges.
# 2. [2013. 10. 10] Added dpx weighted calculation routine

NEUTRALS = ["ALA", "GLY", "ILE", "LEU", "MET", "PRO", "VAL",
			"ASN", "CYS", "GLN", "SER", "THR", "PHE", "TRP",
			"TYR", "ARG", "HIS"]
POSITIVES = ["ARG", "LYS"]
NEGATIVES = ["ASP", "GLU"]

def get_res_list_by_charge(opts = {})
	# Option parsing
	opts = {:threshold => 10.0, :rsa_file => nil, :dpx => false }.merge(opts)
	protein = opts.fetch(:protein)
	charge = opts.fetch(:charge)
	rsa_file = opts.fetch(:rsa_file)
	threshold = opts.fetch(:threshold)
	
	rsa_file = nil if opts.fetch(:dpx) == true
	
	list = Array.new
	neutrals = NEUTRALS
	positives = POSITIVES
	negatives = NEGATIVES
	
  pos = protein.find_residue{|r| r.resName == "LYS" || r.resName == "ARG"}
  neg = protein.find_residue{|r| r.resName == "ASP" || r.resName == "GLU"}
  
  if (rsa_file != nil)
    rsa = parse_rsa_file rsa_file
    del = Array.new
    rsa.each do |l|
      del << l[0] if l[1] < threshold
    end
    pos.delete_if {|res| del.include?(res.resSeq)}
    neg.delete_if {|res| del.include?(res.resSeq)}
  end
  
  case charge.upcase
  when "P" then return pos
  when "N" then return neg
  else return []
  end
  
  # residues = protein.residues.sort_by {|x| x.resSeq}
  # if (rsa_file != nil)
  #   rsa = parse_rsa_file rsa_file
  #   del = Array.new
  #   rsa.each do |l|
  #     del << l[0] if l[1] < threshold
  #   end
  #   residues.delete_if {|res| del.include?(res.resSeq)}
  # end
  # 
  # residues.each do |res|
  #     res_name = res.resName.upcase
  #   case charge.upcase
  #     when "P" then list << res if (res_name == "ARG") || (res_name == "LYS")
  #     when "N" then list << res if (res_name == "ASP") || (res_name == "GLU")
  #     else list << res if neutrals.include?(res_name)
  #   end
  # end
  # return list # returns an array of residues
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

def get_dpx_weighted_vector_sum protein, dpx_file, res_list
  dpx = parse_dpx_protein dpx_file, res_list
  vectors = res_list.map {|res| get_normalized_residue_vector(protein, res)}
  output = Array.new
  vectors.each do |vector|
    i = vectors.index(vector)
    output[i] = 9 * vector / (dpx[i] * dpx[i])
  end
  output.inject {|sum, x| sum + x}
end

def parse_dpx_protein pdb_file, res_list
  protein = Bio::PDB.new(IO.read(pdb_file))
  output = Array.new
  res_list.each do |res|
    resseq, resname = res.resSeq, res.resName
    atoms = protein.find_atom {|a| a.resSeq == resseq && a.resName == resname}
    value = 0
    atoms.each do |atom|
      value = value + atom.tempFactor
    end
    output << value / atoms.size
  end
  return output
end

def parse_rsa_file rsa_file
	# FIXED 2013-07-24
	# as rsa file is a width-fixed text file, data should be parsed as line[0..2]
	# rather than line.split(" ")
	values = Array.new
	IO.foreach(rsa_file) do |line|
		if line[0..2] == "RES"
			# FIXED 2013-08-12
			# as some rsa files are delimited by space(s), data should not be
			# parsed using text position.
			
			#id = line[10..12].strip.to_i
			#rsa = line[23..27].strip.to_i
			#chain_id = line[8]
			data = line[9..-1].split(" ").map{|x| x.to_f}
			id, rsa = data[0], data[2]
			values << [id.to_i, rsa]
		end
	end
	a = values.sort_by {|x| x[0]}
	return a
end

def calc_procd_score(opts = {})
	# option parsing
	opts = {:threshold => 10.0, :rsa_file=>nil, :mutations=>nil, :dpx=>false }.merge(opts) # default threshold is 10.0
	pdb_file = opts.fetch(:pdb_file) # name of pdb file 
	rsa_file = opts.fetch(:rsa_file) # name of rsa file
	threshold = opts.fetch(:threshold) # RSA threshold
	mutations = opts.fetch(:mutations) # mutations
	dpx = opts.fetch(:dpx)
	if dpx
	  rsa_file = nil
	  mutations = nil
	  dpx_file = pdb_file[0..-5] + "_dpx.pdb"
	end
	
	# parse protein
	begin
		protein = Bio::PDB.new(IO.read(pdb_file))
	rescue
		return {:error_code => 1}
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
	
	if pos_size == 0
		pp = 0
	else
	  if dpx
	    pp = get_dpx_weighted_vector_sum(protein, dpx_file, pos)
	  else
		  pp = get_vector_sum(protein, pos)
		end
	end
	
	if neg_size == 0
		nn = 0
	else
	  if dpx
	    nn = get_dpx_weighted_vector_sum(protein, dpx_file, neg)
	  else
		  nn = get_vector_sum(protein, neg)
		end
	end

	# mutation
	if mutations != nil
    case mutations.upcase
    when "SUC"
      lys = protein.find_residue{|r| r.resName == "LYS"}
      lys_vec = get_vector_sum(protein, lys)
      lys_num = lys.size
      pp = pp - lys_vec
      nn = nn + lys_vec
      pos_size = pos_size - lys_num
      neg_size = neg_size + lys_num
    when "AMI"
      pp = pp + nn
      nn = 0
      pos_size = pos_size + neg_size
      neg_size = 0
    else
		  mutations.split(" ").each do |mutation|
			  # parse mutation, mutation should have this form, E421Q
			  orig = Bio::AminoAcid::Data::NAMES[mutation[0]].upcase # 3-letter code
			  mut = Bio::AminoAcid::Data::NAMES[mutation[-1]].upcase # 3-letter code
			  site = mutation[1..-2].to_i # mutation site number
		
			  res = protein.find_residue{|res| res.resName == orig && res.resSeq == site}[0]
			  if res == nil
			    return {:error_code => 2}
			  end
		
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
    end
	end
	
	if (pp == 0) || (nn == 0)
		angle = 0
		score = 0
	else
		angle = calc_angle pp, nn
		score = pp.magnitude/nn.magnitude
	end
	
	# Result output
	# puts "Positive Residues: #{pos_size}"
	# puts "Negative Residues: #{neg_size}"
	# puts "Magnitude of P sum vector: #{pp.magnitude.round(3).to_s}"
	# puts "Magnitude of N sum vector: #{nn.magnitude.round(3).to_s}"
	# puts "Angle between two vectors: #{angle.round(2).to_s}"
	# puts "ProCD Score: #{(pp.magnitude/nn.magnitude).round(3)}"
	
	values = {:pos_size => pos_size, :neg_size => neg_size,
			:pp_length => pp.magnitude.round(3),
			:nn_length => nn.magnitude.round(3),
			:angle => angle.round(2),
			:score => score.round(3),
			:error_code => 0}
	return values
end

def calc_angle a, b
	Math::acos((a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (a.magnitude * b.magnitude)) * 180 / Math::PI
end

def auto_mutate(opts={})
	# Option parsing
	opts = {:threshold => 10.0, :rsa_file=>nil, :chain_id=>nil}.merge(opts)
	pdb_file = opts.fetch(:pdb_file) # name of pdb file 
	chain_id = opts.fetch(:chain_id) # chain ID
	rsa_file = opts.fetch(:rsa_file) # name of rsa file
	threshold = opts.fetch(:threshold).to_f # RSA threshold

	rsa_option = case threshold
		when nil
			"none"
		else
			"filter"
	end

	# calculate original procd score
	orig = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file, :threshold => threshold)
	
	time_stamp = Time.new.strftime("%Y%m%d_%H%M%S_%L")
	output_file_name = "#{pdb_file}_#{time_stamp}.txt"
	begin
		if chain_id == nil
			protein = Bio::PDB.new(IO.read(pdb_file))
		else
			protein = Bio::PDB.new(IO.read(pdb_file)).models[0][chain_id]
		end
	rescue
		puts "Protein file error!"
		return nil
	end
	residues = protein.residues
	
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
		res_name = res.resName.upcase
		res_name_cap = res.resName.capitalize
		res_name_1 = Bio::AminoAcid.to_1(res_name_cap) # 1-letter amino acid code
				
		pos_mutation = "#{res_name_1}#{id.to_s}K"
		neg_mutation = "#{res_name_1}#{id.to_s}D"
		del_mutation = "#{res_name_1}#{id.to_s}A"
		case res_name
			when "ASP", "GLU"
				del = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
					:threshold => threshold, :mutations=>del_mutation)
				pos = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
					:threshold => threshold, :mutations=>pos_mutation)
				output << ["[#{res_name_cap}]#{id.to_s}:#{chain_id}", "delete", "#{del.fetch(:pos_size)}",
				"#{del.fetch(:neg_size)}", "#{del.fetch(:pp_length)}", "#{del.fetch(:nn_length)}", "#{del.fetch(:score)}"]
				output << ["[#{res_name_cap}]#{id.to_s}:#{chain_id}", "neg->pos", "#{pos.fetch(:pos_size)}",
				"#{pos.fetch(:neg_size)}", "#{pos.fetch(:pp_length)}", "#{pos.fetch(:nn_length)}", "#{pos.fetch(:score)}"]
			when "LYS", "ARG"
				del = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
					:threshold => threshold, :mutations=>del_mutation)
				neg = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
					:threshold => threshold, :mutations=>neg_mutation)
				output << ["[#{res_name_cap}]#{id.to_s}:#{chain_id}", "delete", "#{del.fetch(:pos_size)}",
				"#{del.fetch(:neg_size)}", "#{del.fetch(:pp_length)}", "#{del.fetch(:nn_length)}", "#{del.fetch(:score)}"]
				output << ["[#{res_name_cap}]#{id.to_s}:#{chain_id}", "pos->neg", "#{neg.fetch(:pos_size)}",
				"#{neg.fetch(:neg_size)}", "#{neg.fetch(:pp_length)}", "#{neg.fetch(:nn_length)}", "#{neg.fetch(:score)}"]
			else
				pos = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
					:threshold => threshold, :mutations=>pos_mutation)
				neg = calc_procd_score(:pdb_file => pdb_file, :rsa_file => rsa_file,
					:threshold => threshold, :mutations=>neg_mutation)
				output << ["[#{res_name_cap}]#{id.to_s}:#{chain_id}", "positive", "#{pos.fetch(:pos_size)}",
				"#{pos.fetch(:neg_size)}", "#{pos.fetch(:pp_length)}", "#{pos.fetch(:nn_length)}", "#{pos.fetch(:score)}"]
				output << ["[#{res_name_cap}]#{id.to_s}:#{chain_id}", "negative", "#{neg.fetch(:pos_size)}",
				"#{neg.fetch(:neg_size)}", "#{neg.fetch(:pp_length)}", "#{neg.fetch(:nn_length)}", "#{neg.fetch(:score)}"]
		end
	end
	
	CSV.open(output_file_name, "wb") do |csv|
		output.each do |line|
			csv << line
		end
	end
end