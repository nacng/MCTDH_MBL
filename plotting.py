from spin_glass_inset_weak import q_inset_weak
from spin_glass_inset_intermediate import q_inset_int
from mutual_info_inset_unscaled import mi_inset_unscaled
from dynamics_shortTimes import short_time_wi
from dynamics_scaled import scaled_wi
from dynamics_unscaled import fixed_coupling
from dynamics_strong import dynamics_strong


if __name__=='__main__':
    q_inset_weak()
    q_inset_int()
    mi_inset_unscaled()
    short_time_wi()
    dynamics_strong()
    scaled_wi()
    fixed_coupling()
