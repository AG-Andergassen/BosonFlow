#pragma once

#include <models/chain_hubbard.h>
#include <models/square_hubbard.h>
#include <models/fcc_hubbard.h>
#include <models/triangular_hubbard.h>
#include <models/square_hubbard_holstein.h>
#include <models/hubbard_atom.h>
#include <models/anderson_impurity.h>
#include <models/square_hubbard_long_range.h>
#include <models/anderson_impurity_holstein.h>
#include <models/square_hubbard_peierls.h>

#include <params_physical.h>

#define WITH_MODELS(F) F(THE_MODEL)
