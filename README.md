# Distributed-Feedback-Fiber-Laser
This Matlab code computes the characteristics of a rare-earth doped distributed feedback (DFB) fibre laser
Distributed feedback Fibre laser is a type of short cavity fibre laser. As the name says it, its 
feedback is distributed throughout the cavity by a pi phase shifted fibre Bragg grating printed in 
the core of the Erbium-Ytterbium co-doped gain medium.
The main advantage of a DFB fibre laser is its single frequency operation, which is very important for
applications in telecommunication and sensing.
The two counter-propagating fields inside the DFB fibre laser periodic structure is represented by 
two coupled first order differential equations
These equations represent a 2 point boundary value problem and are solve with a shooting algorithm
implemented in Matlab to yield characteristics of the DFB fibre laser.
Using this Matlab model, output power, internal gain, threshold and other characteristics of the laser
can be obtained.
