.. libtsm documentation master file, created by
   sphinx-quickstart on Thu May 27 14:48:45 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

libtsm
==================================
``libtsm`` is a Python toolbox for time-scale modification (TSM) and pitch-shifting.

Details and example application:

https://www.audiolabs-erlangen.de/resources/MIR/2021-DAFX-AdaptivePitchShifting

The toolbox is based on a re-implementation of the
`Matlab TSM toolbox <https://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/>`_ by Jonathan Driedger and Meinard Müller.

If you use this toolbox, please consider the following references:

.. [#] Sebastian Rosenzweig, Simon Schwär, Jonathan Driedger, and Meinard Müller: Adaptive Pitch-Shifting with Applications to Intonation Adjustment in A Cappella Recordings, Proceedings of the International Conference on Digital Audio Effects (DAFx), 2021.

.. [#] Jonathan Driedger and Meinard Müller: TSM Toolbox: MATLAB Implementations of Time-Scale Modification Algorithms. In Proceedings of the International Conference on Digital Audio Effects (DAFx): 249–256, 2014.

.. [#] Jonathan Driedger and Meinard Müller: A Review on Time-Scale Modification of Music Signals. Applied Sciences, 6(2): 57–82, 2016.

.. [#] Jonathan Driedger, Meinard Müller, and Sebastian Ewert: Improving Time-Scale Modification of Music Signals using Harmonic-Percussive Separation. IEEE Signal Processing Letters, 21(1): 105–109, 2014.


.. toctree::
    :hidden:

    getting_started



.. toctree::
    :caption: API Documentation
    :maxdepth: 1
    :hidden:

    index_tsm
    index_pitchshift
    index_utils
