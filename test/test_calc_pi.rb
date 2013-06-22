require 'test/unit'
require './lib/calc_pi'

class TestCalcPI < Test::Unit::TestCase

	def test_count_res
		seq = get_seq_from_pdb("./test/2D9Q.pdb", "A")
		assert_equal count_res(seq, "Asp"), 4
		assert_equal count_res(seq, "CYS"), 0
	end
  
  def test_percentage_res
    seq = get_seq_from_pdb("./test/2D9Q.pdb", "A")
    assert_equal percentage_res(seq).fetch("Ala"), 10.92
    assert_block do
      percentage_res(seq).values.any? {|num| num >= 0}
    end
  end

	def test_calc_nq
		assert_equal calc_nq(6.5, 4, 9, 5, 3, 5, 4, 5).round(2), -2.79
	end

	def test_calc_isoelectric_point
		assert_equal calc_isoelectric_point("./test/2D9Q.pdb", "A"), 5.60
	end
end