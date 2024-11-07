module mod_periodic_boundary
    use constants
    implicit none
    
        
    contains
    function relative_vector(L,v,u)
        real(wp)::L,v(3),u(3),relative_vector(3)
        integer(i4b)::k
        relative_vector=u-v
        where (relative_vector>= L/2)
            relative_vector=relative_vector-L
        else where (relative_vector< -L/2)
            relative_vector=relative_vector+L
        end where
    end function 
end module mod_periodic_boundary