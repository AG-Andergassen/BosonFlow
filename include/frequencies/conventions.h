#include <mymath.h>


/*
The four frequency convention here is this one (used in the IR)
                            +--------+
                 -w2   -->--|        |-->--   w3
                            |   V    |
                  w1   --<--|        |--<--  -w4
                            +--------+
The signs were chosen in the IR side of things so that energy
conservation is only a sum w1 + w2 + w3 + w4 = 0

Everywhere else in the code, this corresponds to 
   w1_in = -w2
   w1_out = w3
   w2_in = -w4
   w2_out = w1

Remark: in the code, the arguments correspond to the integers of 
the frequencies, which sometimes give the extraneously looking 1s.
*/

namespace Symmetrised4P
{
    void pp_to_full_freqs(int W, int w, int wp, int &w1, int &w2, int &w3, int &w4)
    {
	const int div2_ceil_W = div2_ceil(W);
	const int div2_floor_W = W - div2_ceil_W;
    
	w1 = div2_floor_W - wp - 1;
	w2 = -w - div2_ceil_W;
	w3 = wp + div2_ceil_W;
	w4 = -div2_floor_W + w + 1;
    }

    void full_freqs_to_pp(int w2, int w3, int w4, int &W, int &w, int &wp)
    {
	W = -w2 - w4 + 1;
	w = -w2 - div2_ceil(W);
	wp = w3 - div2_ceil(W);
    }

    void ph_to_full_freqs(int W, int w, int wp, int &w1, int &w2, int &w3, int &w4)
    {
	const int div2_ceil_W = div2_ceil(W);
	const int div2_floor_W = W - div2_ceil_W;

	w1 = -div2_floor_W + wp;
	w2 = -w + div2_floor_W;
	w3 = div2_ceil_W + w;
	w4 = -wp - div2_ceil_W;
    }

    void full_freqs_to_ph(int w2, int w3, int w4, int &W, int &w, int &wp)
    {
	const int div2_ceil_W = div2_ceil(W);
	const int div2_floor_W = W - div2_ceil_W;
    
	W = w3 + w2;
	w = -w2 + div2_floor_W;
	wp = -w4 - div2_ceil_W;   
    }

    void xph_to_full_freqs(int W, int w, int wp, int &w1, int &w2, int &w3, int &w4)
    {
	const int div2_ceil_W = div2_ceil(W);
	const int div2_floor_W = W - div2_ceil_W;

	w1 = div2_ceil_W + w;
	w2 = -w + div2_floor_W;
	w3 = -div2_floor_W + wp;
	w4 = -wp - div2_ceil_W;
    }

    void full_freqs_to_xph(int w2, int w3, int w4, int &W, int &w, int &wp)
    {
	const int div2_ceil_W = div2_ceil(W);
	const int div2_floor_W = W - div2_ceil_W;
    
	W = -w4 - w3;
	w = -w2 + div2_floor_W; // w_pietro = -w2 + ceil - floor + floor = w_ours - floor + ceil
	wp = -w4 - div2_ceil_W; // wp_pietro = -w4 - floor - ceil + ceil = wp_ours + ceil - floor 
    }
}


namespace Symmetrised3P
{
    /*
                            +------
                 -w2   -->--|      \   
                            |lambda ~~>~~ W =  -w2 - w4
                 -w4   -->--|      /
                            +------
       Energy conservation: - (-w2) + - (-w4) + W = w2 + w4 + W = 0
       Remark: the symbols below are the matsubara integers, rather than the frequency itself.
       So this would mean that the arguments sum up -1
    */
    void pp_to_full_freqs(int W, int w, int &w2, int &w4)
    {
	// W = W
	const int div2_ceil_W = div2_ceil(W);

	w2 = - w - div2_ceil_W;
	w4 = - W + w + div2_ceil_W - 1;
	// note: w2 + w4 + W = 1
    }

    void full_freqs_to_pp(int W, int w2, int &w)
    {
	w = -w2 - div2_ceil(W);
    }

    /*
                            +------
                 -w2   -->--|      \   
                            |lambda ~~>~~ W = w3 + w2
                  w3   --<--|      /
                            +------
       - (-w2) + w3 + W = w2 + w3 + w4 = 0
       Remark: the symbols below are the matsubara integers, rather than the frequency itself.
       They don't sum up to one in this case however, due to cancelation from -w2 and w3
    */
    void ph_to_full_freqs(int W, int w, int &w2, int &w3)
    {
	const int div2_floor_W = div2_floor(W);

	w2 = -w + div2_floor_W;
	w3 = w - div2_floor_W + W;//erroneous "-W" replaced by "+W" here
    }

    void full_freqs_to_ph(int W, int w2, int &w)
    {
	w = -w2 + div2_floor(W);
    }


    /*
                            +------
                 -w2   -->--|      \   
                            |lambda ~~>~~ W = w1 + w2
                  w1   --<--|      /
                            +------
       w1 + w2 + W = 0
       everything same as ph, but with w3 swapped with w1. It's redundant, but we keep
       the functions below for clarity
    */
    void xph_to_full_freqs(int W, int w, int &w1, int &w2)
    {
	const int div2_floor_W = div2_floor(W);
	w2 = -w + div2_floor_W;
	w1 = w - div2_floor_W + W;//erroneous "-W" replaced by "+W" here
    }

    void full_freqs_to_xph(int w2, int W, int &w)
    {
	w = -w2 + div2_floor(W);
    }
}
