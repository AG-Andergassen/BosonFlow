#include <frg/observables_common.h>
#include <frg/sbe/state.h>

// a global container that contains all observables to be tracked and outputted
template <typename Model, typename state_t >
class observables_frg_sbe_t : public observables_frg_common_t<Model, state_t >
{
 public:    
 observables_frg_sbe_t():observables_frg_common_t<Model, state_t >() {};

    void update(const state_t &state, double t, bool is_final = false) override;    

    static void SetObservablesToTrack();

    static void SetToTrackAllAvailableObservables();
};

// include implementation
#include <../src/frg/sbe/observables.tpp>
