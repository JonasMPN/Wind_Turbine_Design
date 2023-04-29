function FAST_object = blade(FAST_object, aero_data, structure_data, rotor_precone)
    hub_radius = FAST_object.Nacelle.Housing.Diameter/2;
    n_positions = size(aero_data.BlSpn);
    
    %% change Radius
    FAST_object.Blade.Radius = aero_data.BlSpn+hub_radius;
    
    %% change Chord
    FAST_object.Blade.Chord = aero_data.BlChord;
    
    %% change Twist and necessary pitch control values
    fine_pitch = aero_data.BlTwist(end);
    FAST_object.Blade.Twist = aero_data.BlTwist-fine_pitch;
    FAST_object.Control.Pitch.Fine = fine_pitch;
    FAST_object.Control.Pitch.Min = -2+fine_pitch;
    disp("The twist distribution has been reset to be zero at the tip. " + ...
        "The fine pitch was adjusted accordingly.")
    
    %% change NFoil
    FAST_object.Blade.NFoil = double(1:n_positions)';
    
    %% change IFoil
    FAST_object.Blade.IFoil = double(1:n_positions)';
    
    %% change Mass density
    FAST_object.Blade.Mass = structure_data.BMassDen;
    
    %% change flap stiffness
    FAST_object.Blade.EIflap = structure_data.FlpStff;
    
    %% change edge stiffness
    FAST_object.Blade.EIedge = structure_data.EdgStff;
    
    %% change pitch axis
    FAST_object.Blade.PitchAxis = structure_data.PitchAxis;
    
    %% change thickness
    FAST_object.Blade.Thickness = ones(n_positions(1),1);

    %% change rotor precone
    FAST_object.Blade.Cone = rotor_precone;
end

