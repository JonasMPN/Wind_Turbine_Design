function FAST_object = nacelle(FAST_object, nacelle_data)
    FAST_object.Nacelle.Hub.ShaftTilt = nacelle_data.shaft_tilt;
    FAST_object.Nacelle.Hub.Mass = nacelle_data.hub_mass;
    FAST_object.Nacelle.Housing.Diameter = 2*nacelle_data.hub_radius;
    FAST_object.Nacelle.Housing.Mass = nacelle_data.nacelle_mass;
end

