# edited at 04OCT13
# edited at 06OCT13
# edited at 07OCT13

require 'bio'
require 'bio/db/pdb'
require 'FileUtils'

# READ PDB FILE 
def read_pdb(pdb_file)
  text = IO.read(pdb_file).split("\n").delete_if {|x| x =~ /^HETATM/}.join("\n")
  pdb = Bio::PDB.new(text)

  if pdb.record("EXPDTA").to_s.include?("SOLUTION NMR") == true 
    pdb = pdb.models[0]
  else
  pdb = pdb
  end
  return pdb
end

# ROTATION
def center_of_mass(pdb_file)
  com = read_pdb(pdb_file).centreOfGravity
  return com
end

def move_pdb_to_origin(pdb_file)
  pdb = read_pdb(pdb_file)
  com = center_of_mass(pdb_file)

  pdb.atoms.each do |line|
    line.x = (line.x - com.x).round(3)
    line.y = (line.y - com.y).round(3)
    line.z = (line.z - com.z).round(3)
  end
  return pdb
end

def rot_x(pdb, x)
  xr = x * Math::PI / 180

  pdb.atoms.each do |line|
  # r = radius / t = theta / p = phi
    if (line.y < 0 and line.z < 0) or (line.y < 0 and line.z > 0)
      r = Math.sqrt(line.x ** 2 + line.y ** 2 + line.z ** 2)
      t = Math.acos(line.x / r)
      p = Math.atan(line.z / line.y)
      
      line.x = (r * Math.cos(t)).round(3)
      line.y = (r * Math.sin(t) * Math.cos(Math::PI + p + xr)).round(3)
      line.z = (r * Math.sin(t) * Math.sin(Math::PI + p + xr)).round(3)
    else
      r = Math.sqrt(line.x ** 2 + line.y ** 2 + line.z ** 2)
      t = Math.acos(line.x / r)
      p = Math.atan(line.z / line.y)
      
      line.x = (r * Math.cos(t)).round(3)
      line.y = (r * Math.sin(t) * Math.cos(p + xr)).round(3)
      line.z = (r * Math.sin(t) * Math.sin(p + xr)).round(3)
    end
  end
  return pdb
end

def rot_y(pdb, y)
  yr = y * Math::PI / 180

  pdb.atoms.each do |line|
  # r = radius / t = theta / p = phi
    if (line.z < 0 and line.x < 0) or (line.z < 0 and line.x > 0)
      r = Math.sqrt(line.x ** 2 + line.y ** 2 + line.z ** 2)
      t = Math.acos(line.y / r)
      p = Math.atan(line.x / line.z)
      
      line.x = (r * Math.sin(t) * Math.sin(Math::PI + p + yr)).round(3)
      line.y = (r * Math.cos(t)).round(3)
      line.z = (r * Math.sin(t) * Math.cos(Math::PI + p + yr)).round(3)
    else
      r = Math.sqrt(line.x ** 2 + line.y ** 2 + line.z ** 2)
      t = Math.acos(line.y / r)
      p = Math.atan(line.x / line.z)
      
      line.x = (r * Math.sin(t) * Math.sin(p + yr)).round(3)
      line.y = (r * Math.cos(t)).round(3)
      line.z = (r * Math.sin(t) * Math.cos(p + yr)).round(3)
    end
  end
  return pdb
end

def rot_z(pdb, z)
  zr = z * Math::PI / 180

  pdb.atoms.each do |line|
  # r = radius / t = theta / p = phi
    if (line.x < 0 and line.y < 0) or (line.x < 0 and line.y > 0)
      r = Math.sqrt(line.x ** 2 + line.y ** 2 + line.z ** 2)
      t = Math.acos(line.z / r)
      p = Math.atan(line.y / line.x)
      
      line.x = (r * Math.sin(t) * Math.cos(Math::PI + p + zr)).round(3)
      line.y = (r * Math.sin(t) * Math.sin(Math::PI + p + zr)).round(3)
      line.z = (r * Math.cos(t)).round(3)
    else
      r = Math.sqrt(line.x ** 2 + line.y ** 2 + line.z ** 2)
      t = Math.acos(line.z / r)
      p = Math.atan(line.y / line.x)
      
      line.x = (r * Math.sin(t) * Math.cos(p + zr)).round(3)
      line.y = (r * Math.sin(t) * Math.sin(p + zr)).round(3)
      line.z = (r * Math.cos(t)).round(3)
    end
  end
  return pdb
end

# def xyz(pdb_file)
#   time = Time.new.strftime("%Y%m%d-%H%M%S-%L") 
#   x = rand(-360..360)
#   y = rand(-360..360)
#   z = rand(-360..360)
#   xyz = z(y(x(pdb_file, x), y), z)
#   File.open("#{time}_#{pdb_file[0..3].upcase}.pdb_opt_x=#{x},y=#{y},z=#{z}.pdb", "wb") do |file|
#     xyz.atoms.each do |line|
#       file << line
#     end
#   end
#   puts "Input PDB file: #{pdb_file}"
#   puts "Rotation Angle on X-axis: #{x} degrees"
#   puts "Rotation Angle on Y-axis: #{y} degrees"
#   puts "Rotation Angle on Z-axis: #{z} degrees"
#   puts ""
#   return "#{time}_#{pdb_file[0..3].upcase}.pdb_opt_x=#{x},y=#{y},z=#{z}.pdb"
# end

# GROMACS PROCESS
def solvate_protein(protein, output_name)
  if Dir.exist?("temp")
    FileUtils.rm_rf Dir.glob("temp/*")
  else
    Dir.mkdir("temp")
  end
  
  Dir.mkdir("temp_solvated_pdb") if Dir.exist?("temp_solvated_pdb") == false

  File.open("temp/temp.pdb", "wb") do |file|
    protein.atoms.each {|line| file << line}
  end
  
  if File.exist?("temp/temp.pdb")
    `pdb2gmx -f temp/temp.pdb -o temp/conf.gro -p temp/topol.top -i temp/porse.itp -ff amber96 -water tip3p -ignh -missing > /dev/null`
  else
    return "Cannot find .pdb file"
  end
  
  if File.exist?("temp/conf.gro")
    `editconf -f temp/conf.gro -o temp/box.gro -d 1.0 > /dev/null`
  else
    return "Cannot find .gro file"
  end
  
  if File.exist?("temp/box.gro")
    `genbox -cp temp/box.gro -cs -p temp/topol.top -o temp_solvated_pdb/#{output_name}.pdb > /dev/null`
  else
    return "Cannot find .gro file"
  end
  puts "#{Time.new} => gromacs process finished."
end

def make_25_proteins(protein)
  FileUtils.rm_rf Dir.glob("temp_solvated_pdb/*")
  solvate_protein(protein, "output0")
  
  (1..24).each do |i|
    output_name = "output#{i.to_s}"
    case i % 3
    when 0 
      p = rot_x(protein, (i/3)*40)
    when 1
      p = rot_y(protein, (i+2)/3*40)
    when 2
      p = rot_z(protein, (i+1)/3*40)
    end
    solvate_protein(p, output_name)
  end
  puts "#{Time.new} => make 25 proteins finished."
end

# DPX CALCULATOR
def parse_gromacs_pdb(gro_pdb)
  pdb = read_pdb(gro_pdb)
  atom = Array.new
  water = Array.new

  pdb.atoms.each do |el|
    name, resName = el[1], el[3]

    if name !~ /^H.*/ and name !~ /^\dH.*/
       unless resName == "SOL"
         atom << el
       else
         water << el
       end
    else
       nil
    end
  end
  coordination = {"atom" => atom, "water" => water}
  return coordination
end

def get_distance(a, b)
  dist = Math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2).round(3)
  return dist
end

def remove_non_bulk_water(gro_pdb)
  atom = parse_gromacs_pdb(gro_pdb)["atom"]
  water = parse_gromacs_pdb(gro_pdb)["water"]
  
  residue = Hash.new{|k,v| k[v] = Array.new}
  (0..atom.size - 1).each do |i|
    if atom[i].resSeq != atom[i + 1]
      residue["#{atom[i].resName}:#{atom[i].resSeq}"] << [atom[i].x, atom[i].y, atom[i].z]
    else
      nil
    end
  end
  residue_center = Hash.new
  residue.map{|k,v| residue_center[k] = [(v.map{|e| e[0]}.reduce(:+) / v.size).round(3), (v.map{|e| e[1]}.reduce(:+) / v.size).round(3), (v.map{|e| e[2]}.reduce(:+) / v.size).round(3)]}

  new_water = Array.new
  residue_center.values.each do |rc|
    water.each do |e|
      if get_distance(rc, e.to_a) < 8
        new_water << e
      else
        nil
      end
    end
  end

  new_water = new_water.uniq
    dist = Hash.new{|k,v| k[v] = Array.new}
  new_water.each do |el|
    new_water.each do |ele|
      dist[el] << get_distance(el.to_a,ele.to_a) if get_distance(el.to_a,ele.to_a) != 0.0
    end
  end

  close_dist = Hash.new
  dist.map{|k,v| close_dist[k] = v.count{|x| x < 4.2}}
  close_dist.each do |k,v|
    if v < 3
      new_water = new_water.delete_if{|e| e == k}
    else
      nil
    end
  end
  return new_water.uniq
end

def min_dist_atom_to_water(gro_pdb)
  coordination_of_atom = parse_gromacs_pdb(gro_pdb)["atom"]
  coordination_of_water = remove_non_bulk_water(gro_pdb)

  dist = Hash.new{|k,v| k[v] = Array.new}
  coordination_of_atom.each do |el|
    el.name = "OXT" if el.name == "OC1"
    el.name = "O" if el.name == "OC2"
    el.name = "CD1" if el.name == "CD"
    coordination_of_water.each do |ele|
      dist[[el.serial, el.name, el.resName, el.resSeq]] << get_distance(el.to_a, ele.to_a)
    end
  end

  min_dist = Hash.new
  dist.map{|k,v| min_dist[k] = v.min}
  return min_dist
end

def calc_avg_dpx
  pdbs = Dir["temp_solvated_pdb/output*.pdb"]
  dpx = min_dist_atom_to_water(pdbs[0])
  keys = dpx.keys
  print "#"
  (1..24).each do |i|
    vals = min_dist_atom_to_water(pdbs[i])
    dpx.keys.each do |k|
      dpx[k] = dpx[k] + vals[k]
    end
  print "#"
  end
  puts "\ndpx calculation is finished!"

  keys.each {|k| dpx[k] = (dpx[k]/25).round(3)}
  return dpx
end

def make_dpx_pdb_file(pdb_file)
  current_time = Time.new.strftime("%Y%m%d-%H%M%S-%L") 
  atom_dist = calc_avg_dpx
  pdb = read_pdb(pdb_file)
  heavy_atom = Array.new
  hydrogen = Array.new
  water = Array.new
  
  # prepare output file name
  output_name = File.dirname(pdb_file) + "/" + File.basename(pdb_file, File.extname(pdb_file)) + "_dpx.pdb"

  pdb.atoms.each do |el|
    name, resName = el[1], el[3]
    name = "CD1" if name == "CD"

  unless resName == "SOL"
      if name !~ /^H.*/ and name !~ /^\dH.*/
        heavy_atom << el
      else
        hydrogen << el
      end
    else
      water << el
    end
  end

  atom_dist.each do |k,v|
    heavy_atom.each do |ele|
      ele.name = "CD1" if ele.name == "CD"
      
      if ele.resSeq == k[3] && ele.name == k[1]
        ele.tempFactor = v
      else
        nil
      end
    end
  end

  File.open(output_name, "wb") do |file|
    heavy_atom.each do |line|
      file << line
    end
  end
  
  # remove temporary files
  FileUtils.rm_r(%w("temp_solvated_pdb" "temp"), :force => true)
end

def dpx_run(pdb_file)
  start_time = Time.new
  puts "#{start_time} => dpx calculation started."
  # prepare 25 solvated proteins
  make_25_proteins(move_pdb_to_origin(pdb_file))
  
  # Using 25 files, calculate dpx and write the result into new pdb file
  make_dpx_pdb_file(pdb_file)
  
  end_time = Time.new
  
  puts ""
  puts "#{start_time} => dpx calculation started."
  puts "#{end_time} => dpx calculation finished."
  puts ""
  puts "total running time => #{(end_time - start_time)}"
  puts ""
  puts "#{pdb_file}_dpx.pdb successfully created!"
end