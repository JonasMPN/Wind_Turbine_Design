function FAST_object = nacelle(FAST_object, shaft_tilt, hub_radius, hub_mass, nacelle_mass)
    FAST_object.Nacelle.Hub.ShaftTilt = shaft_tilt;
    FAST_object.Nacelle.Hub.Mass = hub_mass;
    FAST_object.Nacelle.Housing.Diameter = 2*hub_radius;
    FAST_object.Nacelle.Housing.Mass = nacelle_mass;
end

