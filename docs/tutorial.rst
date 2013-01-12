IMUSim Tutorial
===============

Martin Ling, May 2011

The aim of this tutorial is to quickly introduce some basic use of `IMUSim
<http://www.imusim.org/>`_. The software can do much more than what is
covered here, but working through the basics will introduce you to the
structure of the software and how to use it. Once you have worked through
this tutorial, the `API documentation <http://www.imusim.org/docs/api/>`_
should hopefully make a bit more sense to you, and you can refer to that
for further details.

If you haven't installed IMUSim yet, see the :doc:`Installation instructions
<installation>`.

Introduction
------------

IMUSim is written in the `Python <http;//www.python.org/>`_ programming
language, and the normal way to make use of it is by writing Python. If you
don't know Python, don't worry! It's a very simple language and you can
probably pick it up as you go along. If you'd like to learn some Python before
you start, see the `Python tutorial <http://docs.python.org/tutorial/>`_.

Python is an interpreted language, so you can use it interactively. This is
very useful for experimenting and exploring. You can use the standard Python
interpreter to work with IMUSim. However we highly recommend using the
`IPython shell <http://ipython.scipy.org/>`_, which provides an enhanced
interactive environment with tab completion, improved command history, easy
documentation lookup, and other features.

Once you have installed IMUSim, start IPython from the command line with::

    $ ipython -pylab -wthread

The -pylab option imports the Pylab environment, which makes available lots
of additional functions and classes for scientific computing. If you've used
MATLAB before, you'll find most of the same facilities. The -wthread option
starts a wxWidgets event loop, which will make the 3D visualisations used
by IMUSim work correctly.

Now you can follow along with the commands in the tutorial. Each line
starting with '>>>' is a command to enter into the interpreter.

Getting started
---------------

To use the functions and classes provided by IMUSim, you need to import them
into your Python namespace. IMUSim provides an imusim.all package which makes
it easy to import everything from IMUSim at once::

    >>> from imusim.all import *

This might take a few seconds to run since IMUSim loads a lot of additional
libraries.

Now that all the functions and classes from the IMUSim API are available to us,
we can create a new IMUSim simulation::

    >>> sim = Simulation()

This constructs an instance of the :api:`Simulation` class, and names it
`sim`.

Now we need something to simulate. Let's start by simulating a single IMU.
Before we can simulate an IMU, we need to define the trajectory it will follow
through time, since this determines the values its sensors will measure.
IMUSim provides various ways to define trajectories, which we will explore
later. To start with, let's create a random trajectory::

    >>> trajectory = RandomTrajectory()

Now we'll create a simulated IMU. There are several IMU models. We'll start
with an idealised one, :api:`IdealIMU`. To create a simulated IMU, we need to
specify a simulation and trajectory::

    >>> imu = IdealIMU(sim, trajectory)

The IMU includes simulated sensors and other components, but we haven't told
it what to do with them yet. We need to define the behaviour of the IMU, which
simulates the code which would be running on it. Typically this would include
sampling the sensors and running a filter or some other processing on each set
of samples as they are generated. Fortunately, we have a standard behaviour
that does just that, called :api:`BasicIMUBehaviour`. It takes several
optional arguments which let you define the processing to run, but the only
essential ones are the IMU to run on, and the interval at which to sample::

    >>> samplingPeriod = 0.01
    >>> behaviour = BasicIMUBehaviour(imu, samplingPeriod)

Now our IMU will sample all its sensors every 0.01 seconds, i.e. at 100Hz. We
haven't specified any processing for it to run on the data, but all the
sampled sensor values will be recorded.

We're almost ready to run our simulation. The only thing we haven't set is
the start time. A trajectory is only defined for a given time period, so
we can only simulate an IMU within that period. We'll set the time within
the simulation to the first time at which our trajectory is defined::

    >>> sim.time = trajectory.startTime

Now we can run our simulation by calling its :api:`Simulation.run` method,
which takes as an argument the time at which to stop simulating. We'll run
the simulation up to the last time at which the trajectory is defined::

    >>> sim.run(trajectory.endTime)
    Simulating...
    Simulated 0.1s of 1.8s (  5%). Estimated time remaining 0.4s
    ...
    Simulated 1.7s of 1.8s ( 95%). Estimated time remaining 0.0s
    Simulation complete.
    Simulated 1.8 seconds in 0.4 seconds.

Now we can look at the results. Let's plot the accelerometer samples from
the IMU::

    >>> plot(imu.accelerometer.rawMeasurements)

Of course, we should label our graph appropriately::

    >>> title("Accelerometer Readings")
    >>> xlabel("Time (s)")
    >>> ylabel("Acceleration (m/s^2)")
    >>> legend()

The result is shown below:

.. image:: acceleration-plot.png

Graphs are plotted using the matplotlib library, so see the `matplotlib
documentation <http://matplotlib.sourceforge.net/>`_ for much more information
on plotting options. Note that IMUSim has overridden the normal plot command
to work with its data types. You can still use IMUSim's :api:`plot` in all the
normal ways.

Data types
----------

What did we just pass to :api:`plot`? Let's have a look at it::

    >>> imu.accelerometer.rawMeasurements
    <imusim.utilities.time_series.TimeSeries object at 0xa5358d0>

It's a :api:`TimeSeries` object, one of the basic data types of IMUSim, used
to represent time series data. It has two main attributes, timestamps and
values. The timestamps attribute is an array of time values in ascending
order::

    >>> imu.accelerometer.rawMeasurements.timestamps
    array([ 0.01,  0.02, ...,  1.79,  1.8 ])

These are times at which the samples were taken. The sample values themselves
are found in the values attribute::

    >>> imu.accelerometer.rawMeasurements.values
    array([[  66.705814  ,   55.14694052, ..., -211.04480933, -204.6486176 ],
           [ -93.40026896,  -98.81340505, ..., -161.00547201, -155.16993659],
           [ 116.56420017,   94.92124762, ...,  108.21077098,  117.56964057]])

Note the shape of this array, which is 3xN where N is the number of
timestamps. This is the format in which IMUSim works with arrays of vector
data, indexed first by axis and then by sample number. A single vector would
be represented as a 3x1 array. IMUSim has a :api:`vector` function to
construct these::

    >>> vector(1,2,3)
    array([[ 1.],
           [ 2.],
           [ 3.]])

The array type comes from `NumPy <http://numpy.scipy.org/>`_. If you're not
familiar with NumPy arrays then it may be worth reading the first few sections
of the `NumPy Tutorial <http://www.scipy.org/Tentative_NumPy_Tutorial>`_.

The other important data type is the quaternion, which is a mathematical
construct with four components that IMUSim uses to represent a rotation in 3D
space. IMUSim provides its own :api:`Quaternion` class, which can be
constructed with given values, or created from other rotation representations
such as Euler angle sequences::

    >>> Quaternion(1,0,0,0)
    Quaternion(1.0, 0.0, 0.0, 0.0)
    >>> Quaternion.fromEuler((45, 10, 30), order='zyx')
    Quaternion(0.89763565965718295, 0.20599112279858978, 0.17644656798009617, 0.34739673068127142)

They can also be converted to other representations such as rotation
matrices::

    >>> Quaternion(1,0,0,0).toMatrix()
    matrix([[ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.]])

If you are not familiar with the use of quaternions to represent rotations,
you can see `this page
<http://www.genesis3d.com/~kdtop/Quaternions-UsingToRepresentRotation.htm>`_
for an brief mathematical introduction, or for an in-depth treatment consult
the book "Quaternions and Rotation Sequences" by Jack B. Kuipers. However,
there is no need to understand quaternion mathematics to use the simulator,
you just need to know that each quaternion represents a rotation.

We introduced a time series of vectors above, but the :api:`TimeSeries` class
can also be used with scalar or quaternion values. Let's look at a quaternion
time series: the rotation key frames of the random trajectory we used
earlier::

    >>> trajectory.rotationKeyFrames
    <imusim.utilities.time_series.TimeSeries object at 0xa8e5f10>
    >>> trajectory.rotationKeyFrames.values
    QuaternionArray(array([[-0.04667991, -0.82763194,  0.29852366, -0.47300103],
                           [-0.10730387, -0.8172798 ,  0.33822558, -0.45402981],
                           ..., 
                           [ 0.40666661, -0.04250812,  0.80062722,  0.43796277],
                           [ 0.42667426, -0.01498163,  0.82309247,  0.37449087]]))

Arrays of quaternions are stored using the special :api:`QuaternionArray`
class, which wraps an Nx4 NumPy array of the component values. Quaternion
arrays provide support for applying quaternion math operations over the whole
array.

As with the vector time series, we can plot a quaternion time series simply
by passing it to :api:`plot`, and adding labels as appropriate::

    >>> figure() # Create a new figure
    >>> plot(trajectory.rotationKeyFrames)
    >>> title("Quaternion components")
    >>> xlabel("Time (s)")
    >>> ylabel("Component value")
    >>> legend()

The result is shown below:

.. image:: quaternion-plot.png

Trajectories
------------

Now that we've looked at the vector and quaternion types, we're ready to
understand the trajectory classes. A trajectory defines the path of an object
through space, and also its changing rotation, over time. If that object is
a simulated IMU, we can obtain the (ideal) readings of its sensors -
accelerometer, gyroscope and magnetometer - from the acceleration, rotational
rate and rotation of its trajectory at a given time.

So a trajectory needs to provide position and rotation, and their first and
second derivatives, at any given time. Let's look at the methods used to
obtain these. We'll look at all the trajectory parameters at the starting
time of the trajectory, which it can tell us::

    >>> t = trajectory.startTime
    >>> t
    3.8146809461460811

The position, velocity and acceleration methods of the trajectory provide
vector values at time t::

    >>> trajectory.position(t)
    array([[-10.36337587],
           [  4.63926506],
           [ -0.17801693]])
    >>> trajectory.velocity(t)
    array([[ 30.79525389],
           [-20.9180481 ],
           [  2.68236355]])
    >>> trajectory.acceleration(t)
    array([[ 178.30674569],
           [ -15.11472827],
           [  15.54901256]])

The rotation at time t is a quaternion, but its derivatives - rotational
velocity and acceleration - are vectors::

    >>> trajectory.rotation(t)
    Quaternion(-0.046679914286751743, -0.82763194419093666, 0.29852365607750142, -0.47300103310567193)
    >>> trajectory.rotationalVelocity(t)
    array([[-2.97192064],
           [ 2.97060751],
           [-7.32688967]])
    >>> trajectory.rotationalAcceleration(t)
    array([[ -8.46813312],
           [ 19.43475152],
           [-31.28760834]])

The last thing a trajectory needs to provide is its end time::

    >>> trajectory.endTime
    5.6159546456299427

That's all a trajectory needs to do. Any object which can implement the
methods we've just looked at can be used as a trajectory by IMUSim. Note that
in this case the whole trajectory was defined in advance, but it is possible
to use a trajectory which is defined as a simulation progresses, e.g. by
simulating the effect of some control system. The simulator will only call
the trajectory methods for a time when all events prior to that time have been
simulated.

IMUSim provides a number of trajectory classes to help you define realistic
trajectories easily. See the documentation for the :api:`trajectories` module.

Environment models
------------------

The sensor readings of an IMU don't just depend on its trajectory, they also
depend on the environment. Accelerometers sense gravity, and magnetometers
sense magnetic field, both of which vary with position. We may also want
to simulate radio transmissions from a wireless IMU, the propagation of which
will depend on its surroundings. All of these environmental considerations are
described by an Environment object. Each simulation has an associated
Environment, which in turn has models for each of the environmental aspects
of the simulation::

    >>> sim.environment
    <imusim.environment.base.Environment object at 0x4405910>
    >>> sim.environment.magneticField
    <imusim.environment.magnetic_fields.EarthMagneticField object at 0xa557750>
    >>> sim.environment.gravitationalField
    <imusim.environment.gravity.ConstantGravitationalField object at 0xa111690>
    >>> sim.environment.radioEnvironment
    <imusim.environment.radio_environment.IdealRadioEnvironment object at 0xa6b11d0>

This is the default environment, which was created automatically with our
simulation. The gravitationalField and magneticField attributes are vector
field models, i.e. vector quantities that depend on position and time. We can
evaluate the field values at a given position by calling the models with a
position vector and time::

    >>> p = trajectory.position(t)
    >>> sim.environment.gravitationalField(p, t)
    array([[ 0.  ],
           [ 0.  ],
           [ 9.81]])
    >>> sim.environment.magneticField(p, t)
    array([[  1.71010072e-05],
           [  0.00000000e+00],
           [  4.69846310e-05]])

These are the values of standard Earth gravity, 9.81m/s^2, and the approximate
Earth magnetic field (in Teslas) in Edinburgh, UK where IMUSim was written.
You can get figures for your own location from the `International Geomagnetic
Reference Field model <http://www.ngdc.noaa.gov/geomagmodels/IGRFWMM.jsp>`_
and pass them to the :api:`EarthMagneticField` constructor.

Both the gravity and magnetic field models are instances of the VectorField
class. This has several subclasses for modelling fields in various ways. E.g.
the :api:`SolenoidMagneticField` class models the magnetic field around a
single ideal solenoid. More complex fields can be modelled by superposition of
multiple solenoids. Or, if you have a set of field samples at known positions
you can interpolate a field model from them using the
:api:`NaturalNeighbourInterpolatedField` class. At present, none of the field
models included with the simulator are time-varying, but time varying fields
are supported if you create one.

Sensor and IMU models
---------------------

Of course, the measurements from real sensors depend on much more than just
their trajectory and environment. They suffer from noise, bias, misalignment,
cross-axis sensitivity and many other effects. To acheive a realistic
simulation we need to model these. IMUSim includes generic models for
imperfect sensors with various parameters, and also specific models of some
real sensor components, derived from measurements and datasheet information.

It's not just the sensors themselves which affect the measurements, though.
The voltages produced by the sensors must be converted to digital values using
an analogue to digital converter (ADC), which clips and quantises the values
as well as adding its own noise. Plus, the samples are never taken at the
exact times we wanted, because of the inevitable inaccuracy of the IMU's
hardware timers. Also, multiple sensors are often multiplexed into the same
ADC and sampled in sequence, rather than simultaneously as we would prefer.
So to account for all these problems, we also need to model the ADC and timers
of an IMU. IMUSim includes some generic parametric models for these components.

All of these components can be brought together to create a model of a
specific IMU design. The :api:`IdealIMU` we used earlier is an example, with
ideal models for all its components. IMUSim also includes a model of the real
Orient-3 IMU developed at Edinburgh, :api:`Orient3IMU`. Any IMU model can be
modified easily by simply assigning different components to its relevant
attributes before simulation. Or, you can write your own sensor and/or IMU
classes.

IMU calibration
---------------

Let's start putting together a more realistic simulation. We'll start by
taking an instance of the Orient-3 IMU model::

    >>> imu = Orient3IMU()

Notice that this time we're not passing a simulation and trajectory to the IMU
constructor. We're actually going to put this IMU through several simulations
on different trajectories, so we'll assign these parameters later.

The IMU model we've just created includes many imperfections. Its outputs are
noisy, and it has biases and scaling factors between the real values it senses
and the numbers which come out of its ADC. So just like a real IMU, we need to
calibrate it before it can be used for measurements. We calibrate an IMU by
putting it through a set of controlled trajectories, looking at the sensor
values obtained, and working out how to transform the raw measurements into
calibrated ones with meaningful units. In IMUSim, procedures for doing this are
implemented as :api:`Calibrator` classes. A calibrator takes an :api:`IMU`
instance and calibrates it for a given :api:`Environment`. It returns a set of
:api:`SensorCalibration` objects which transform measurements from each sensor
of the IMU into calibrated values.

The current release includes a single calibrator, :api:`ScaleAndOffsetCalibrator`.
As explained in its documentation, it implements a simple calibration procedure
which can be followed by hand for real IMUs. It simulates testing the device in
18 different scenarios - static with each accelerometer axis pointing up and
down, then with each magnetometer axis pointing north and south, and then
rotating at a known speed around each axis (e.g. on a rotation stage). Based on
the results, it fits a constant scale and offset to each axis of each sensor on
the IMU.

To initialise the calibrator we need to provide several parameters - the
environment to calibrate in, the number of samples to take in each scenario,
the sampling period, and the angular speed in rad/s for the rotating tests::

    >>> env = Environment()
    >>> samples = 1000
    >>> rotationalVelocity = 20
    >>> calibrator = ScaleAndOffsetCalibrator(env, samples, samplingPeriod, rotationalVelocity)

Now we can calibrate our IMU::

    >>> calibration = calibrator.calibrate(imu)

This may take a little while to run - it's running the IMU through 18 different
simulations and then fitting the calibration.

When it's done, the result is a Python dictionary object, mapping from each
:api:`Sensor` object on the IMU to a :api:`ScaleAndOffsetCalibration` object
with the calibration parameters for that sensor. Let's look at the scale and
offset parameters obtained for one of the sensors::

    >>> calibration[imu.accelerometer]
    <imusim.algorithms.calibration.ScaleAndOffsetCalibration object at 0x9fb3490>
    >>> calibration[imu.accelerometer].scale
    array([[ 0.03672461],
           [ 0.03693735],
           [ 0.03260453]])
    >>> calibration[imu.accelerometer].offset
    array([[-115.60203444],
           [  53.67615273],
           [   6.02866548]])

Each is a vector containing values for each axis. You don't need to use these
values yourself, because each calibration object also has an
:api:`SensorCalibration.apply` method that transforms raw measured values into
calibrated ones using these parameters. And you can pass the whole calibration
dictionary as an argument to :api:`BasicIMUBehaviour` too, which will
then create a `calibratedMeasurements` :api:`TimeSeries` attribute on each IMU
sensor and update it as it runs.

There are many different approaches to calibrating IMUs, and the one
implemented by the :api:`ScaleAndOffsetCalibrator` class is just an example.
Others may appear in future releases, and you can implement your own too, using
the current code as an example to guide you.

Of course, we could also cheat by taking the correct transforms directly from
the IMU models. But modelling the calibration procedure accurately, including
its limitations, helps us to acheive realistic results.

Importing motion capture data
-----------------------------

In our first simulation example we used a randomly generated trajectory. These
can be useful for testing, but in most cases it's more useful to use a
realistic trajectory. For some applications it may be possible to define a
suitable trajectory mathematically, in which case a trajectory class can be
written directly. For most applications however, and particularly for complex
movements such as those of a human, the best way to acheive realistic
simulations is to use trajectories derived from motion capture data.

There are two main types of motion capture data formats. Some give positions,
and perhaps also rotations, of independent markers. Others, especially those
used to store human motion capture data, store the movements of a connected
skeletal model. IMUSim supports both types.

Marker-based captures are represented by the :api:`MarkerCapture` class.
The :api:`imusim.io.qualisys_tsv` and :api:`imusim.io.vicon_csv` modules
provide loaders which will produce :api:`MarkerCapture` objects from data
exported by Qualisys and Vicon capture systems.

Skeletal capture data is represented by the :api:`SampledBodyModel` class.
The :api:`imusim.io.bvh` module provides importers and exporters for skeletal
capture data in the common BVH format, and the :api:`imusim.io.asf_amc` module
supports loading ASF/AMC file combinations.

If you have capture data in another format, it is probably relatively
straightforward to write a loader for it, or to write or find a program which
will translate it into one of the supported formats.

Let's load a skeletal capture file. We'll use the `walk.bvh` file which can be
downloaded from `here <http://www.imusim.org/docs/examples/walk.bvh>`_ (97KB)::

    >>> model = loadBVHFile('walk.bvh', CM_TO_M_CONVERSION)

Note the conversion factor, which is simply a multiplier for all the distance
measurements in the file::

    >>> CM_TO_M_CONVERSION
    0.01

This file has measurements in centimetres, but we want metres, so we multiply
by 0.01. The BVH format doesn't define the units of measurement, and there are
files in circulation with various different units. IMUSim uses SI units
throughout, so ensure you convert as necessary.

3D visualisation
----------------

At this point it would be useful to see what's going on in the capture file
we've just imported. IMUSim includes some basic 3D visualisation facilities to
help with this, in the :api:`imusim.visualisation.rendering` module, based on
`MayaVi <http://code.enthought.com/projects/mayavi/>`_.

Visualisations are implemented by renderers, which implement the interface
defined by the abstract :api:`AnimatedRenderer` class. Some example renderers
are included in IMUSim, including the :api:`BodyModelRenderer` which renders
skeletal body models like the one we just imported::

    >>> renderer = BodyModelRenderer(model)

There are two ways to use renderers - either creating a video file using the
:api:`createVideo` function, or creating an interactive animation using the
:api:`InteractiveAnimation` class. We'll create an interactive animation to
display the movements of our skeletal model::

    >>> start = model.startTime
    >>> end = model.endTime
    >>> animation = InteractiveAnimation(start, end, renderer)

You should now have two new windows - one showing a visualisation of the lower
body of a human, and one for controlling the animation, as seen below:

.. image:: walk-visualisation.png

The view can be rotated, moved and zoomed using the mouse. For more
information on the interaction controls see the `Mayavi documentation
<http://github.enthought.com/mayavi/mayavi/application.html#interaction-with-the-scene>`_.
The animation can be started and stopped using the control window. As you will
see, the BVH file we have imported is a capture of the lower body of a walking
human subject.

Rigid body systems
------------------

The skeletal model we have loaded above is an example of a jointed rigid body
model. The model is a tree of joints, each of which has a fixed location in
the local co-ordinate frame of its parent. There may also be points in the
model which are not joints, but lie in a joint co-ordinate frame - these 
represent e.g. the ends of the outermost segments. These two roles are
represented by the :api:`Point` and :api:`Joint` classes. Pure :api:`Point`
and :api:`Joint` objects do not have any movement information, and can be used
to represent just the structure of a model.

The :api:`SampledBodyModel` class is a subclass of :api:`Joint`, representing
a root joint (one with no parent) whose position and rotations are known in
the global co-ordinate frame at discrete times. As a :api:`Joint`, it has a
name, and a set of children, which are also joints, having their own children::

    >>> model.name
    'root'
    >>> model.children
    [<imusim.trajectories.rigid_body.SampledJoint object at 0xa842450>,
     <imusim.trajectories.rigid_body.SampledJoint object at 0x4612a50>]
    >>> [child.name for child in model.children]
    ['rfemur', 'lfemur']

It is also a subclass of :api:`SampledPositionTrajectory` and
:api:`SampledRotationTrajectory`, the trajectory classes for representing
sampled position and rotation data. We can ask it for its position and
rotation at a given time::

    >>> t = model.startTime
    >>> model.position(t)
    array([[-1.51953],
           [-0.06939],
           [-0.97422]])
    >>> model.rotation(t)
    Quaternion(0.50494478394432252, -0.53571366129112219, 0.45998422065524641, -0.49644350637476148)

However, because this is a sampled trajectory it moves in discrete steps,
jumping from one value to the next at the time of each sample. 

Also, we cannot obtain derivatives::

    >>> model.acceleration(t)
    NotImplementedError: Derivative not available from sampled trajectory.
    Create a splined trajectory to obtain derivatives.

As such, we cannot use this type of trajectory to simulate an IMU.

Trajectories from motion capture data
-------------------------------------

IMUSim supports fitting spline functions to sampled position and rotation data
to obtain continuous, differentiable trajectories. There are several wrapper
classes provided which take sampled capture data and return spline-interpolated
versions. For our body model data we can use the :api:`SplinedBodyModel`
class::

    >>> splinedModel = SplinedBodyModel(model)

To see the difference, let's create another animation. Close the previous
animation windows, and then create a new animation with semi-transparent
red and blue rendererings for the sampled and splined body models::

    >>> sampledRenderer = BodyModelRenderer(model, opacity=0.5, color=(1,0,0))
    >>> splinedRenderer = BodyModelRenderer(splinedModel, opacity=0.5, color=(0,0,1))
    >>> animation = InteractiveAnimation(start, end, sampledRenderer, splinedRenderer)

Slow the animation down to the minimum speed before starting it, and zoom
in for a closer look. Notice how the red sampled model moves in discrete
steps, whilst the blue splined model follows the same path but is smooth and
continuous, as shown below:

.. image:: sampled-vs-splined.gif

Similarly, spline interpolation can be applied to capture data based on
independent, uncoupled markers, by using the :api:`SplinedMarkerCapture`
wrapper on :api:`MarkerCapture` data.

A more realistic simulation
---------------------------

Now we have all the components required for a fairly realistic simulation
of an IMU undergoing a real human motion. We'll create a new simulation using
the environment we calibrated in earlier, assign the Orient-3 IMU we created
and calibrated to it, and attach the IMU to the left foot of our
spline-interpolated walking subject::

    >>> sim = Simulation(environment=env)
    >>> imu.simulation = sim
    >>> imu.trajectory = splinedModel.getJoint('rfoot')

We'll start the simulation from the start time of our splined model's
movements. Note that this will be slightly after those of our sampled model,
because some boundary samples are required at the start and end to
unambiguously fit the spline parameters::

    >>> sim.time = splinedModel.startTime

As before, we'll run a :api:`BasicIMUBehaviour` on the IMU, but this time
we will provide it with the calibration information we obtained earlier, and
also initially set the IMU's local clock correctly::

    >>> BasicIMUBehaviour(imu, samplingPeriod, calibration, initialTime=sim.time)

Now we'll run our simulation to the end of our splined model's movements::

    >>> sim.run(splinedModel.endTime)
    Simulating...
    Simulated 0.2s of 3.9s (  5%). Estimated time remaining 1.3s
    ...
    Simulated 3.9s of 3.9s (100%). Estimated time remaining 0.1s
    Simulation complete.
    Simulated 3.9 seconds in 1.4 seconds.

And plot the calibrated accelerometer readings::

    >>> figure()
    >>> plot(imu.accelerometer.calibratedMeasurements)
    >>> title("Accelerometer Readings")
    >>> xlabel("Time (s)")
    >>> ylabel("Acceleration (m/s^2)")
    >>> legend()

The result is shown below. The magnetometer and gyroscope readings can be
examined in the same way.

.. image:: walking-acceleration-plot.png

Conclusion
----------

In this tutorial we have covered basic use of the IMUSim software, up to the
point of being able to obtain realistic sensor readings. Our last simulation
was based on:

    - a real human motion, imported from motion capture data.
    - an empirically obtained model of a real IMU design, including noise and
      other imperfections.
    - a simulation of a realistic calibration procedure.

And it required just 16 lines of code::

    >>> from imusim.all import *
    >>> samplingPeriod = 0.01
    >>> imu = Orient3IMU()
    >>> env = Environment()
    >>> samples = 1000
    >>> rotationalVelocity = 20
    >>> calibrator = ScaleAndOffsetCalibrator(env, samples, samplingPeriod, rotationalVelocity)
    >>> calibration = calibrator.calibrate(imu)
    >>> model = loadBVHFile('walk.bvh', CM_TO_M_CONVERSION)
    >>> splinedModel = SplinedBodyModel(model)
    >>> sim = Simulation(environment=env)
    >>> imu.simulation = sim
    >>> imu.trajectory = splinedModel.getJoint('rfoot')
    >>> sim.time = splinedModel.startTime
    >>> BasicIMUBehaviour(imu, samplingPeriod, calibration, initialTime=sim.time)
    >>> sim.run(splinedModel.endTime)

Obtaining realistic sensor data in simulations is one of IMUSim's key goals,
but there is also much more it can do. So far, all we have done with our sensor
data is to plot it. IMUSim also includes implementations of many algorithms for
processing IMU sensor data to reconstruct the orientation and position of an
IMU, or the movements of a body model with multiple attached IMUs. Many of
these may be useful as they are in many applications, but IMUSim was also
intended to be used as a platform for developing new algorithms and techniques,
and provides many mathematical primitives and useful features to help with
this.

Now that you've worked through this tutorial, you can refer to the `API
documentation <http://www.imusim.org/docs/api/>`_ documentation for further
details on all the classes and functions provided by IMUSim. If you have
questions, you can ask them on the IMUSim mailing list. Details are on the
`IMUSim home page <http://www.imusim.org/>`_.

We hope you enjoy using IMUSim and find it useful!
