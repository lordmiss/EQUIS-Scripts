$:.unshift File.join(File.dirname(__FILE__), "..", "lib")
require 'tm_calculator'

require 'test/unit'

class TestNAcid < Test::Unit::TestCase
	
	def setup
		@seq = NAcid.new("ggctggtgcaagtcacagacttggctg")
		@empty_seq = NAcid.new("")
	end

	def test_calc_tm
		assert_equal(73.2, @seq.calc_tm[:tm])
		assert_equal(64.29, @seq.calc_tm[:b_tm])
		assert_equal(@seq.calc_tm[:b_tm], @seq.calc_tm(:conc_na=>0.10)[:b_tm])
		assert_equal(78.2, @seq.calc_tm(:conc_na=>0.10)[:tm])
		assert_equal(3, @seq.calc_tm.size)
		assert_equal(1, @empty_seq.calc_tm[:error_code])
	end

	def test_get_f_half
		assert_equal(13, @seq.get_f_half.size)
		assert_equal("g", @seq.get_f_half[0].to_s)
	end

	def test_get_t_half
		assert_equal("acagacttggctg", @seq.get_t_half.to_s)
	end

	def test_is_symmetric?
		assert(!(@seq.is_symmetric?))
	end

	def test_get_symmetric_part
		assert_kind_of(Array, @seq.get_symmetric_part)
		assert_equal(2, @seq.get_symmetric_part.size)
		assert_equal(["caagt", 8], @seq.get_symmetric_part)
	end
	
	def test_find_seq_at_tm
		assert_nil(@seq.find_seq_at_tm(100))
		assert_nil(@seq.find_seq_at_tm(12))
		assert_equal("ggctggtgcaag", @seq.find_seq_at_tm(40))
		assert_equal(20, @seq.find_seq_at_tm(60).size)
	end

	end # of TestNA