module UtilTest

using Test
import KLIEPInference: trimap, itrimap

@testset "trimap" begin
    @test trimap(2, 1) == 1
	@test trimap(1, 2) == 1
	@test trimap(4, 3) == 6

	@test itrimap(1) == (2, 1)
	@test itrimap(6) == (4, 3)
end

end
