=================
A note on units 
=================

The model's numerical work is done in a nanometer, gram, second unit-space to, well, match the scale of the system. This means that we use conversion factors to get from mks units to our nmgs system. Our unit of force is piconewtons (:math:`1 g*nm^2/s^2`). Conversion factors for common units are:

* Joules (:math:`1 kg\,m^2/s^2`) multiply by :math:`10^{21}` to get it in :math:`pN\,nm`
* Poise (:math:`1 g/cm\,s`) multiply by :math:`10^{-7}` to get it in :math:`g/nm\,s`
