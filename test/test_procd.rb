require 'test/unit'
require './lib/procd'

class TestProCD < Test::Unit::TestCase
  @@protein = Bio::PDB.new(IO.read("test/2D9Q_A.pdb"))
  @@residues = @@protein.residues
  @@rsa_file = "test/2D9Q_A.rsa"
  def test_get_res_list_by_charge
    # without rsa filter
    assert_equal get_res_list_by_charge(:protein=>@@protein, :charge=>"P").size, 9
    assert_equal get_res_list_by_charge(:protein=>@@protein, :charge=>"N").size, 13
    # with rsa filter
    assert_equal get_res_list_by_charge(:protein=>@@protein, :rsa_file=>@@rsa_file, :threshold=>10.0, :charge=>"N").size, 12
    # wrong option (other than p and n)
    assert_equal get_res_list_by_charge(:protein=>@@protein, :charge=>"a").size, 0
  end
  
  def test_get_normalized_residue_vector
    # 1. result with any residue should be 1 as it is normalized
    assert_block do
      @@residues.any? {|res| get_normalized_residue_vector(@@protein, res).norm == 1}
    end
  end
  
  def test_get_vector_sum
    assert_equal get_vector_sum(@@protein, get_res_list_by_charge(:protein=>@@protein, :charge=>"P")).norm.round(3), 4.311
  end
  
  def test_parse_rsa_file
    # 1. 2D9Q has 168 residues
    assert_equal parse_rsa_file(@@rsa_file).size, 168
    # 2. all element of returned array should have exactly 2 elements
    assert_block do
      parse_rsa_file(@@rsa_file).any? {|e| e.size == 2}
    end
  end
  
  def test_calc_procd_score
    # 1. basic calculation
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdb", :rsa_file => @@rsa_file,
	  :threshold=>10.0).fetch(:score), 1.039
    # 2. calculation with 1-point mutation
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdb", :rsa_file => @@rsa_file,
	  :threshold=>10.0, :mutations => "D109A").fetch(:score), 1.304
    # 3. protein file error case (error code 1)
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdx", :rsa_file => @@rsa_file,
	  :threshold=>10.0, :mutations => "D109A").fetch(:error_code), 1
    # 4. mutation input error (error code 2)
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdb", :rsa_file => @@rsa_file,
	  :threshold=>10.0, :mutations => "Q109A").fetch(:error_code), 2
    # 5. when threshold value is big, the result should be zero
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdb", :rsa_file => @@rsa_file,
	  :threshold=>100.0).fetch(:score), 0
    # 6. mutations option1 : suc (succinylation)
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdb", :rsa_file => @@rsa_file,
	  :threshold=>10.0, :mutations => "suc").fetch(:score), 0.296
    assert_equal calc_procd_score(:pdb_file => "test/2D9Q_A.pdb", :rsa_file => @@rsa_file,
	  :threshold=>10.0, :mutations => "suc").fetch(:pos_size), 5
  end
end