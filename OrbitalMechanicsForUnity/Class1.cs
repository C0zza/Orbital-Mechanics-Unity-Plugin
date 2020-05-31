using System;
using UnityEngine;
using UnityEditor;

namespace OrbitalMechanicsForUnity  // Plugin to implement orbital mechanics into the Unity 3D Engine
{
    [Serializable]
    public class Body : MonoBehaviour   
    {
        public double   bodyDiameter = 0;
        public double   mass = 5972000000000000000000000d;    // Mass in kg of the body
        public double   hillSphereRadius;                     // The radius of this influence
        public bool     DrawHillSphere = true;                  // Whether to draw the body's area of dominant gravitational influence on other objects.
        public Body     Primary = null;                         // The body being orbited (if any)

        [SerializeField]
        Orbit       orbit = new Orbit();  // Orbital elements
        DVector3    position;          // Universe position
        bool        accelerating = false;
        DVector3    velocityVector;
        public void UpdateBody()    // Update body parameters/ details
        {
            if (UniverseManager.GetTimeScale() != 0 && Primary) // If time has not stopped and is orbiting something
            {   
                if(accelerating)
                {
                    if(velocityVector == null) // if velocity vector has been initialised
                    {
                        velocityVector = OrbitalMechanics.VelocityVector(orbit, this);
                        DVector3.LogDVector3(velocityVector);
                    }

                    velocityVector += new DVector3(gameObject.transform.forward * 10) / mass; // Apply thrust
                                                                                                                          // Calculate centripetal force (gravity of primary)
                    DVector3 centripetalForce = OrbitalMechanics.PositionVector(this).Normalized() * -Math.Pow((OrbitalMechanics.Velocity(Primary.mass, Primary.position,
                                                                                    position, orbit.SemiMajorAxis)), 2d) / DVector3.Distance(position, Primary.position);
                    DVector3.LogDVector3(centripetalForce);

                    velocityVector += centripetalForce; // Apply centripetal force
                    velocityVector *= UniverseManager.GetTimeScale();   // Adjust to current time scale (high time scales not recommended for accelerating)

                    position += velocityVector * Time.deltaTime;  // Update position with new velocity vector
                    DVector3.LogDVector3(position);

                    orbit = OrbitalMechanics.OrbitalElementsFromStateVectors(OrbitalMechanics.PositionVector(this), velocityVector, Primary.mass, this);    // Update orbital elements
                }
                else
                {
                    // Update position in orbit by passing time passed this frame
                    orbit.TrueAnomaly = OrbitalMechanics.TrueAnomalyAfterTime(orbit.SemiMajorAxis, orbit.Eccentricity, orbit.TrueAnomaly,
                                                                   Time.deltaTime * (float)UniverseManager.GetTimeScale(), Primary.mass);

                    if (orbit.TrueAnomaly > 360) orbit.TrueAnomaly -= 360;      // Clamp the true anomaly 0-360 to prevent the value from
                    else if (orbit.TrueAnomaly < 0) orbit.TrueAnomaly += 360;   // getting too large, potentially causing data issues

                    position = GetWorldPosition();  // Update universe position
                }
            }                              
                
        }

        private void Awake()
        {
            GameObject bodyObject = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            bodyObject.transform.localScale = new Vector3((float)bodyDiameter, (float)bodyDiameter, (float)bodyDiameter);
            bodyObject.transform.parent = gameObject.transform;
            bodyObject.transform.localPosition = new Vector3(0, 0, 0);
        }

        private void Start()    // Standard Unity start function
        {
            position = GetWorldPosition();  // Initialise universe position

            //if (UniverseManager.scaled && Primary)  // If universe is currently being scaled and this body is orbiting something
            //{                                                             // Set game position to new universe position
            //    gameObject.transform.localPosition = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit)) /
            //                      (float)OrbitalMechanics.AstronomicalUnit * (float)UniverseManager.GetUniverseScale();
            //}

            gameObject.transform.position = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit));
            //if (Primary) gameObject.transform.position += Primary.transform.position;
            gameObject.transform.position -= DVector3.GetFloatVector(UniverseManager.GetCurrentOrigin());
        }
        private void Update()   // Standard Unity update function
        {
            UpdateBody();

            //if (UniverseManager.scaled && Primary)  // If universe is currently being scaled and this body is orbiting something
            //{                                                             // Set game position to new universe position
            //    gameObject.transform.localPosition = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit)) /
            //                      (float)OrbitalMechanics.AstronomicalUnit * (float)UniverseManager.GetUniverseScale();
            //}

            gameObject.transform.position = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit));
            //if (Primary) gameObject.transform.position += Primary.transform.position;
            gameObject.transform.position -= DVector3.GetFloatVector(UniverseManager.GetCurrentOrigin());
        }
        public void StartAccelerating()
        {
            accelerating = true;
        }
        public void StopAccelerating()
        {
            accelerating = false;
            velocityVector = null;
        }
        public void AddOrbitingBody(double mass = 1)
        {
            GameObject newBody = Instantiate(new GameObject(), transform);
            newBody.transform.parent = transform;
            newBody.name = "Body";
            GameObject tempSphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            tempSphere.transform.parent = newBody.transform;
            tempSphere.transform.localPosition = new Vector3(0, 0, 0);
            tempSphere.transform.localScale = new Vector3(0.1f, 0.1f, 0.1f);

            Body newBodyComponent = newBody.AddComponent<Body>();
            if(Primary)
            {
                newBodyComponent.Orbit().SemiMajorAxis = (hillSphereRadius / UniverseManager.GetUniverseScale()) / 2d;
                Debug.Log(hillSphereRadius);
                newBodyComponent.mass = 73476730900000000000000d;
            }
            else
            {
                newBodyComponent.mass = 597200000000000000000000d;
                newBodyComponent.Orbit().SemiMajorAxis = OrbitalMechanics.AstronomicalUnit;
            }

            newBodyComponent.Primary = this;
        }
        public DVector3 GetPosition()   // Get universe position of object
        {
            return position;
        }
        public DVector3 GetWorldPosition()  // Calculate and get universe position
        {
            if(Primary)     // If body is orbiting something.
            {
                if(Primary.GetPosition() == null)   // If primary's position member is not initialised
                {
                    return OrbitalMechanics.PositionInOrbit(orbit) + OrbitalMechanics.PositionInOrbit(Primary.Orbit()); // Calculate both positions in
                }                                                                                                       // orbit and return their sum
                else
                {
                    return OrbitalMechanics.PositionInOrbit(orbit) + Primary.GetPosition(); // Calculate this body's position and sum with primary position
                }
            }
            else // Body does not orbit anything
            {
                return new DVector3(transform.position) * OrbitalMechanics.AstronomicalUnit; // Convert and return game position to universe position
            }
        }
        public Orbit Orbit()    // Get this body's orbit/ orbital elements
        {
            return orbit;
        }
        private void OnDrawGizmos()
        {
            if(DrawHillSphere)
            {
                Gizmos.color = Color.blue;
                Gizmos.DrawWireSphere(gameObject.transform.position, (float)(hillSphereRadius / OrbitalMechanics.AstronomicalUnit));
            }
        }
    }       // Main component of the plugin, provides the behaviours of a satellite in space

    [CustomEditor(typeof(Body))]    // informs Unity which component it should act as an editor for
    [CanEditMultipleObjects]        // Tells Unity that multiple objects can be selected and edited at the same time
    public class BodyEditor : Editor
    {
        SerializedProperty mass, orbit, drawHillSphere, primary, diameter;    // Parameters to display within the inspector
        void OnEnable() // When the body is created
        {
            mass = serializedObject.FindProperty("mass");
            orbit = serializedObject.FindProperty("orbit");
            drawHillSphere = serializedObject.FindProperty("DrawHillSphere");
            primary = serializedObject.FindProperty("Primary");
            diameter = serializedObject.FindProperty("bodyDiameter");
        }
        public override void OnInspectorGUI()   // Update the inspector display when the body is modified
        {
            serializedObject.Update();

            EditorGUILayout.ObjectField(primary);
            EditorGUILayout.PropertyField(mass);
            EditorGUILayout.PropertyField(diameter);

            Body myBody = (Body)target;

            if(myBody.Primary)
            {
                EditorGUILayout.PropertyField(drawHillSphere);
                EditorGUILayout.PropertyField(orbit);
            }

            if(myBody.mass > 1000000)
            {
                if(GUILayout.Button("Add Orbiting Body"))
                {
                    myBody.AddOrbitingBody();
                }
                //if(GUILayout.Button("Focus on Body"))
                //{
                //    UniverseManager.Focus(myBody);
                //}
            }

            myBody.Orbit().ValidateOrbit();

            serializedObject.ApplyModifiedProperties();
        }
        public void OnSceneGUI()    // Update the scene view when body is modified
        {
            var body = target as Body;  // Get reference to body component

            if (!body) return;  // If no body returned, exit OnSceneGUI

            if(body.Primary)    // If the body is orbiting something
            {
                Orbit orbit = body.Orbit(); // Get reference to body's orbit member 
                                                                                                                     // Calculate hill sphere radius
                body.hillSphereRadius = OrbitalMechanics.HillSphereRadius(new DVector3(body.transform.position) * OrbitalMechanics.AstronomicalUnit,
                                   new DVector3(body.Primary.transform.position) * OrbitalMechanics.AstronomicalUnit, body.mass, body.Primary.mass);

                if (orbit.RenderOrbit)   // If orbit needs to be rendered
                {                                                                                                                               
                    if (!Application.isPlaying) // If application is not running
                    {                                                                                // Set body's position in scene
                        body.gameObject.transform.localPosition = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit)) /
                                               (float)OrbitalMechanics.AstronomicalUnit * (float)UniverseManager.GetUniverseScale();
                    }                                                                                

                    Handles.color = Color.magenta;  // Set orbit color to magenta
                                                    // Get first position in the orbit
                    Vector3 v1 = OrbitalMechanics.PositionInOrbitMeasuredFromCenter(0d, orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination,
                                                       orbit.LongitudeOfAscendingNode, orbit.ArgumentOfPeriapsis, body.Primary.transform.position);
                    Vector3 v2;

                    for (int i = 1; i < orbit.OrbitResolution; i++)
                    {                                                                                    // Calculate next position in orbit
                        v2 = OrbitalMechanics.PositionInOrbitMeasuredFromCenter(((360d / orbit.OrbitResolution) * i), orbit.SemiMajorAxis, orbit.Eccentricity,
                            orbit.Inclination, orbit.LongitudeOfAscendingNode, orbit.ArgumentOfPeriapsis, body.Primary.transform.position);

                        Handles.DrawLine(v1, v2);   // Draw a line in the scene view from old point to new point

                        v1 = v2; // Set old point as newest point
                    }
                                                                                            // Draw final line to complete orbit drawing
                    Handles.DrawLine(v1, OrbitalMechanics.PositionInOrbitMeasuredFromCenter(0d, orbit.SemiMajorAxis, orbit.Eccentricity,
                        orbit.Inclination, orbit.LongitudeOfAscendingNode, orbit.ArgumentOfPeriapsis, body.Primary.transform.position));
                }
            }
        }
    }   // Custom editor class for body

    [Serializable]
    public class Orbit
    {
        public double SemiMajorAxis; // The radius of the elipses longest axis
        [Range(0f, 1f)]
        public double Eccentricity; // How much the orbit deviates from a circular shape. 0 = circle
        [Range(0f, 180f)]
        public double Inclination;  // Angle of the orbit measured against the eclpitic plane at the ascending node
        [Range(0f, 360f)]
        public double LongitudeOfAscendingNode; // The rotation of the orbit around the ecliptic plane
        [Range(0f, 360f)]
        public double ArgumentOfPeriapsis;  // The rotation of the orbit around the orbital plane
        [Range(0f, 360f)]
        public double TrueAnomaly;  // The current position within the orbit

        public bool RenderOrbit;
        [Range(0f, 200f)]
        public int OrbitResolution = 30;

        public Orbit()
        {
            SemiMajorAxis = OrbitalMechanics.AstronomicalUnit;
            Eccentricity = 0d;
            Inclination = 0d;
            LongitudeOfAscendingNode = 0d;
            ArgumentOfPeriapsis = 0d;
            TrueAnomaly = 0d;
            RenderOrbit = true;
        }
        public void ValidateOrbit()
        {
            if (Inclination < 10e-8)
            {
                LongitudeOfAscendingNode = 0d;
                if(Eccentricity < 10e-8)
                {
                    ArgumentOfPeriapsis = 0d;
                }
            }
            else
            {
                if (Eccentricity < 10e-8)
                {
                    ArgumentOfPeriapsis = 0d;
                }
            }
        }
    }   // Orbital elements container

    public class UniverseManager : MonoBehaviour // Manage global parameters effecting all bodies
    {
        static double TimeScale = 1d;
        static double UniverseScale = 100d; // 1 Astronomical Unit = universe scale unity units
        public static bool scaled = false;  // 
        public GameObject currentFocus = null;
        static DVector3 currentOrigin;
        Body[] allBodies;

        void Awake()
        {
            currentOrigin = OrbitalMechanics.PositionInOrbit(currentFocus.GetComponent<Body>().Orbit());
        }

        void Start()
        {
            allBodies = FindObjectsOfType<Body>();
            Debug.Log(allBodies.Length);
            //gameObject.transform.position = DVector3.GetFloatVector(currentFocus.GetComponent<Body>().GetWorldPosition());
        }

        void Update()
        {
            if(Vector3.Distance(currentFocus.transform.position, new Vector3(0,0,0)) > 600)
            {
                Debug.Log("Updating Origin");
                //for(int i = 0; i < allBodies.Length; ++i)
                //{
                //    allBodies[i].gameObject.transform.position -= currentFocus.transform.position;
                //}
                currentOrigin = OrbitalMechanics.PositionInOrbit(currentFocus.GetComponent<Body>().Orbit());
            }
        }

        public static void SetTimeScale(double timeScale)
        {
            TimeScale = timeScale;
        }

        public static void IncreaseTimeScale(double increaseBy)
        {
            TimeScale += increaseBy;
        }

        public static double GetTimeScale()
        {
            return TimeScale;
        }

        public static void SetUniverseScale(double scale)
        {
            UniverseScale = scale;
        }

        public static double GetUniverseScale()
        {
            return UniverseScale;
        }

        public void Focus(Body focus)
        {
            GameObject manager = GameObject.Find("UniverseManager");

            if (!currentFocus)
            {
                currentFocus = manager;
                manager.transform.position = new Vector3(0, 0, 0);
            }

            Vector3 movement = currentFocus.transform.position - focus.transform.position;
            manager.transform.Translate(movement);
            currentFocus = focus.gameObject;
        }

        public static DVector3 GetCurrentOrigin()
        {
            return currentOrigin;
        }
    }

    public class OrbitalMechanics   // The mathematics behind orbital mechanics
    {
        public const double GravitationalConstant = 0.0000000000667408d;    // Gravitational constant parameter, used to calculate GM
        public const double AstronomicalUnit = 149597870700d; // Scaled astronomical unit

        public static double RadiusAtPointInOrbit(double semiMajorAxis, double eccentricity, double trueAnomaly)    // Radius of the point in an orbit at the given true anomaly
        {                                                                                                           
            return ((semiMajorAxis * (1d - Math.Pow(eccentricity, 2d))) / (1d + (eccentricity * Math.Cos(trueAnomaly))));
        }

        public static double OrbitalPeriodNegligibleMass(double primaryMass, DVector3 primaryPosition, DVector3 bodyPosition)
        {
            double r = DVector3.Distance(primaryPosition, bodyPosition);  // Distance between primary and body
            return Math.Sqrt((4d * Math.Pow(Math.PI, 2d) * Math.Pow(r, 2d)) / (GravitationalConstant * primaryMass));
        }   // Calculate a body's orbital period when it's mass is negligible to the primary mass

        public static double OrbitalPeriod(double primaryMass, double semiMajorAxis, double bodyMass)  // Calulate a body's orbital period when both have large masses, i.e. star and planet
        {
            return Math.Sqrt((4d * Math.Pow(Math.PI, 2d) * Math.Pow(semiMajorAxis, 3d)) / (GravitationalConstant * (primaryMass + bodyMass)));
        }

        public static double AngularVelocity(double primaryMass, DVector3 primaryPosition, DVector3 bodyPosition) // Calculate a bodies angular velocity in radians.
        {
            double r = DVector3.Distance(primaryPosition, bodyPosition);
            return Math.Sqrt((GravitationalConstant * primaryMass) / Math.Pow(r, 3d));
        }

        public static double TrueAnomalyAfterTime(double semiMajorAxis, double eccentricity, double trueAnomaly, float timePassed, double primaryMass)
        {   // Explanations of these calculations found at http://www.braeunig.us/space/orbmech.htm (Position in an Elliptical Orbit, Problems 4.13 and 4.14)
            double v0 = trueAnomaly * Math.PI / 180d;   // Initial true anomaly in radians
            double e0 = (Math.Atan2(Math.Sqrt((1d - eccentricity) / (1d + eccentricity)) * Math.Sin(v0 / 2d), Math.Cos(v0 / 2d))) * 2d;    // Initial eccentric anomaly
            double m0 = e0 - eccentricity * Math.Sin(e0);   // Initial mean anomaly
            double n = Math.Sqrt((GravitationalConstant * primaryMass) / Math.Pow(semiMajorAxis, 3));   // Mean motion
            double m = m0 + n * (timePassed);   // mean anomaly at new position
            // Calculating new eccentric anomaly
            double a1 = m + eccentricity * Math.Sin(2d);    // first answer to iteration and will save oldest iteration
            double a2 = 0; // Stores newest iteration

            while(true) // Iterate through until accurate to 5 decimal places
            {
                a2 = m + eccentricity * Math.Sin(a1);

                if (Math.Round(a1, 5) == Math.Round(a2, 5)) break;  // Answer is accurate enough, exit while loop
                else a1 = a2;
            }

            double v = 2d * Math.Atan(Math.Sqrt((1d + eccentricity) / (1d - eccentricity)) * Math.Tan(a2 / 2d));    // Calculate final true anomaly in radians
            return v * 180d / Math.PI;  // return v converted into degrees
        }   // Calculate the true anomaly after an amount of time passed
                                                                                                                            
        public static Vector3 PositionInOrbitMeasuredFromCenter(double trueAnomaly, double semiMajorAxis, double eccentricity, double inclination,  // Position in orbit at true anomaly measured from center of orbit instead of primary
                                                                double longitudeOfAscendingNode, double argumentOfPeriapsis, Vector3 offset)
        {
            double sma = SemiMinorAxis(semiMajorAxis, eccentricity);
            double majorTimesMinor = semiMajorAxis * sma;
            double coordDenominator = Math.Pow(Math.Tan(trueAnomaly * Math.PI / 180d), 2);

            float x = (float)(majorTimesMinor / (Math.Sqrt(Math.Pow(sma, 2) + Math.Pow(semiMajorAxis, 2) * coordDenominator)));
            float y = (float)(majorTimesMinor / (Math.Sqrt(Math.Pow(semiMajorAxis, 2) + Math.Pow(sma, 2) / coordDenominator)));

            if (trueAnomaly >= 90d && trueAnomaly <= 270d) x = -x;
            if (trueAnomaly >= 180d && trueAnomaly <= 360d) y = -y;

            Quaternion pointRotation = Quaternion.Euler((float)-longitudeOfAscendingNode * Vector3.up);
            pointRotation *= Quaternion.AngleAxis((float)inclination, OrbitVector90(semiMajorAxis, eccentricity));
            pointRotation *= Quaternion.AngleAxis((float)-argumentOfPeriapsis - 90f, Vector3.up);

            Vector3 periapsisPos = (PositionInOrbit(0d, semiMajorAxis, eccentricity,    // Get periapsis position
                      inclination, longitudeOfAscendingNode, argumentOfPeriapsis) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale());
            Vector3 apoapsisPos = (PositionInOrbit(180d, semiMajorAxis, eccentricity,   // Get Apoapsis position
                      inclination, longitudeOfAscendingNode, argumentOfPeriapsis) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale());

            Vector3 returnVec = (pointRotation * new Vector3((float)y, 0, (float)x)) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale();
            returnVec += offset;
            returnVec += (apoapsisPos + periapsisPos) / 2f;

            return returnVec;
        }

        public static Vector3 PositionInOrbitMeasuredFromCenter(Orbit orbit, Vector3 offset)
        {
            double sma = SemiMinorAxis(orbit.SemiMajorAxis, orbit.Eccentricity); // Calculate value of orbit semi major axis.
            double majorTimesMinor = orbit.SemiMajorAxis * sma;
            double coordDenominator = Math.Pow(Math.Tan(orbit.TrueAnomaly * Math.PI / 180d), 2);

            float x = (float)(majorTimesMinor / (Math.Sqrt(Math.Pow(sma, 2) + Math.Pow(orbit.SemiMajorAxis, 2) * coordDenominator)));
            float y = (float)(majorTimesMinor / (Math.Sqrt(Math.Pow(orbit.SemiMajorAxis, 2) + Math.Pow(sma, 2) / coordDenominator)));

            if (orbit.TrueAnomaly >= 90d && orbit.TrueAnomaly <= 270d) x = -x;
            if (orbit.TrueAnomaly >= 180d && orbit.TrueAnomaly <= 360d)y = -y;
            
            Quaternion pointRotation = Quaternion.Euler((float)-orbit.LongitudeOfAscendingNode * Vector3.up);
            pointRotation *= Quaternion.AngleAxis((float)orbit.Inclination, OrbitVector90(orbit.SemiMajorAxis, orbit.Eccentricity));
            pointRotation *= Quaternion.AngleAxis((float)-orbit.ArgumentOfPeriapsis - 90f, Vector3.up);
            
            Vector3 periapsisPos = (PositionInOrbit(0d, orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode,
                                                 orbit.ArgumentOfPeriapsis) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale());
            Vector3 apoapsisPos = (PositionInOrbit(180d, orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode,
                                                  orbit.ArgumentOfPeriapsis) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale());
            
            Vector3 returnVec = (pointRotation * new Vector3((float)y, 0, (float)x)) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale();
            returnVec += offset;
            returnVec += (apoapsisPos + periapsisPos) / 2f;

            return returnVec;
        }   // Position in orbit at true anomaly measured from center of orbit instead of primary. Used to draw smooth orbits

        public static Vector3 PositionInOrbit(double trueAnomaly, double semiMajorAxis, double eccentricity, double inclination, double longitudeOfAscendingNode,
                                              double argumentOfPeriapsis) // returns vector3 of position in orbit at given true anomaly
        {   // CREATE DOUBLE PRECISION QUATERNION CLASS?
            double a = trueAnomaly * Math.PI / 180.0d;
            double d = RadiusAtPointInOrbit(semiMajorAxis, eccentricity, a);
            double x = -Math.Sin(trueAnomaly * Math.PI / 180.0d);
            double z = Math.Cos(trueAnomaly * Math.PI / 180.0d);

            Quaternion pointRotation = Quaternion.Euler((float)-longitudeOfAscendingNode * Vector3.up);
            pointRotation *= Quaternion.AngleAxis((float)inclination, OrbitVector90(semiMajorAxis, eccentricity));
            pointRotation *= Quaternion.AngleAxis((float)-argumentOfPeriapsis - 90f, Vector3.up);

            return pointRotation * new Vector3((float)x, 0, (float)z).normalized * (float)(d);
        }

        public static DVector3 PositionInOrbit(Orbit orbit) // returns DVector3 of position in orbit
        {   // CREATE DOUBLE PRECISION QUATERNION CLASS?
            double a = orbit.TrueAnomaly * Math.PI / 180.0d;
            double d = RadiusAtPointInOrbit(orbit.SemiMajorAxis, orbit.Eccentricity, a);
            double x = -Math.Sin(orbit.TrueAnomaly * Math.PI / 180.0d);
            double z = Math.Cos(orbit.TrueAnomaly * Math.PI / 180.0d);

            Quaternion pointRotation = Quaternion.Euler((float)-orbit.LongitudeOfAscendingNode * Vector3.up);
            pointRotation *= Quaternion.AngleAxis((float)orbit.Inclination, OrbitVector90(orbit.SemiMajorAxis, orbit.Eccentricity));
            pointRotation *= Quaternion.AngleAxis((float)-orbit.ArgumentOfPeriapsis - 90f, Vector3.up);

            return new DVector3(pointRotation * new Vector3((float)x, 0, (float)z).normalized * (float)(d));
        }

        public static Orbit OrbitalElementsFromStateVectors(DVector3 positionVector, DVector3 velocityVector, double primaryMass, Body body)
        {   // https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

            DVector3 angularMomentum = DVector3.CrossProduct(positionVector, velocityVector);
            DVector3 nodeVector = DVector3.CrossProduct(new DVector3(0d, 1d, 0d), angularMomentum); // Vector along orbital and reference plane

            // Variables to be used throughout function
            Orbit orbit = new Orbit();  // Orbit to be returned
            double smallNumber = 1.0e-6;  // any value less than this in magnitude is considered 0
            double velocityMagnitude = velocityVector.Magnitude();
            double positionMagnitude = positionVector.Magnitude();
            double angularMomentumMag = angularMomentum.Magnitude();
            double GM = GravitationalConstant * primaryMass;    // Gtravity force of body being orbited

            // Calculate eccentricity
            DVector3 eccentricityVector = ((Math.Pow(velocityMagnitude, 2d) - GM / positionMagnitude) * positionVector - // Vector pointing from apoapsis to 
                                             DVector3.DotProduct(positionVector, velocityVector) * velocityVector) / GM; // periapsis with magnitude of eccentricity
            orbit.Eccentricity = eccentricityVector.Magnitude();

            //More variables used for calculations
            double eccentrictyVectorMag = eccentricityVector.Magnitude();
            double nodeVectorMag = nodeVector.Magnitude();
            double mechanicalEnergy = Math.Pow(velocityMagnitude, 2d) / 2d - GM / positionMagnitude;
            //double period = 0d;

            // Calculate orbital period and semi major axis
            if(orbit.Eccentricity != 1)
            {
                orbit.SemiMajorAxis = -GM / (2d * mechanicalEnergy);
                //period = orbit.SemiMajorAxis * (1d - Math.Pow(orbit.Eccentricity, 2d));
            }
            else
            {
                //period = Math.Pow(angularMomentumMag, 2d) / GM;
                orbit.SemiMajorAxis = double.PositiveInfinity;
            }

            // Calculate inclination
            double temp1 = angularMomentum.y / angularMomentumMag;
            orbit.Inclination = Math.Acos(temp1);                       // cos(Angle) = dot(a,b) / (mag(a) * mag(b))
            orbit.Inclination *= 180d / Math.PI; // Convert to degrees

            if(orbit.Inclination < smallNumber) // If orbit is not inclined
            {
                orbit.LongitudeOfAscendingNode = 0d;    // longitude of ascending node is undefined so set to 0 by default

                if(orbit.Eccentricity < smallNumber) orbit.ArgumentOfPeriapsis = 0d;  // for circular orbits, set argument of periapsis to 0 / ascending node by default
                else // Orbit is inclined but not circular
                {
                    orbit.ArgumentOfPeriapsis = Math.Acos(eccentricityVector.x / eccentrictyVectorMag);
                    orbit.ArgumentOfPeriapsis *= 180d / Math.PI;
                    if (eccentricityVector.z < 0d) orbit.ArgumentOfPeriapsis = 360d - orbit.ArgumentOfPeriapsis;
                }
            }
            else
            {
                // Calculate longitude of ascending node
                orbit.LongitudeOfAscendingNode = Math.Acos(nodeVector.x / nodeVectorMag);
                orbit.LongitudeOfAscendingNode *= 180d / Math.PI;
                if (nodeVector.z < 0d) orbit.LongitudeOfAscendingNode = 360d - orbit.LongitudeOfAscendingNode;

                Vector3 rotatedEccentricityVector = DVector3.GetFloatVector(eccentricityVector);
                if (orbit.LongitudeOfAscendingNode > smallNumber)
                {
                    Quaternion argpeRotation = Quaternion.Euler((float)orbit.LongitudeOfAscendingNode * Vector3.up);
                    rotatedEccentricityVector = argpeRotation * DVector3.GetFloatVector(eccentricityVector);
                }

                // Calculate argument of periapsis                                                                                                        
                orbit.ArgumentOfPeriapsis = Math.Acos(DVector3.DotProduct(nodeVector, eccentricityVector) / (nodeVectorMag * eccentrictyVectorMag));
                orbit.ArgumentOfPeriapsis *= 180d / Math.PI;
                if (rotatedEccentricityVector.z < 0d) orbit.ArgumentOfPeriapsis = 360d - orbit.ArgumentOfPeriapsis;
                
                if(orbit.Inclination > 90)
                {
                    orbit.ArgumentOfPeriapsis = 360d - orbit.ArgumentOfPeriapsis;
                }
            }

            // Calculate true anomaly
            if (orbit.Eccentricity < smallNumber)
            {
                orbit.ArgumentOfPeriapsis = 0d;
                if (orbit.Inclination < smallNumber)
                {
                    orbit.TrueAnomaly = Math.Acos((positionVector.x / positionMagnitude));
                    Debug.Log("true anom before val: " + orbit.TrueAnomaly);
                    if (velocityVector.x < 0d) orbit.TrueAnomaly = 2d * Math.PI - orbit.TrueAnomaly;
                    orbit.TrueAnomaly *= 180d / Math.PI;
                }
                else
                {
                    orbit.TrueAnomaly = Math.Acos(DVector3.DotProduct(nodeVector, positionVector) / (nodeVectorMag * positionMagnitude));
                    if (DVector3.DotProduct(nodeVector, velocityVector) < 0d)
                    {
                        orbit.TrueAnomaly = 2d * Math.PI - orbit.TrueAnomaly;
                    }
                    orbit.TrueAnomaly *= 180d / Math.PI;

                    if(double.IsNaN(orbit.TrueAnomaly))
                    {
                        Vector3 temp = Quaternion.Euler((float)orbit.LongitudeOfAscendingNode * Vector3.up) * DVector3.GetFloatVector(positionVector);
                        if (temp.x > 0)
                        {
                            orbit.TrueAnomaly = 0d;
                        }
                        else
                        {
                            orbit.TrueAnomaly = 180d;
                        }
                    }

                }
            }
            else
            {
                orbit.TrueAnomaly = Math.Acos(DVector3.DotProduct(eccentricityVector, positionVector) / (eccentrictyVectorMag * positionMagnitude));
                orbit.TrueAnomaly *= 180d / Math.PI;
                double temp = DVector3.DotProduct(positionVector, velocityVector);
                if (temp > 0d) orbit.TrueAnomaly = 360d - orbit.TrueAnomaly;
            }

            return orbit;
        }

        public static DVector3 PositionVector(Body body)
        {
            return body.Primary.GetWorldPosition() - body.GetWorldPosition();
        }

        public static DVector3 VelocityVector(Orbit orbit, Body body)   // Vector with direction of current flight path and magnitude of velocity
        {
            Vector3 primaryToBody = body.gameObject.transform.position - body.Primary.gameObject.transform.position;    // Position vector of body relative to primary

                                                                                            // Cross vector between primary to body position vector and orbital plane up vector
            Vector3 cross = Vector3.Cross(primaryToBody, OrbitVectorUp(orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode));

            double fpa = FlightPathAngle(orbit.SemiMajorAxis, orbit.Eccentricity, body);    // Calculate flight path angle
            fpa *= 180d / Math.PI;  // Convert to degrees
            if (orbit.TrueAnomaly >= 180d) fpa = -fpa;                                                  // Flight path angle rotation applied to be applied to cross vector
            Quaternion rotation = Quaternion.Euler((float)fpa * OrbitVectorUp(orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode));

            Vector3 final = rotation * cross;

            return new DVector3(final.normalized) * Velocity(body.Primary.mass, body.Primary.GetWorldPosition(), body.GetWorldPosition(), orbit.SemiMajorAxis); // Return normalised velocity direction vector multiplied by velocity magnitude
        }

        public static double Velocity(double primaryMass, DVector3 primaryPos, DVector3 bodyPos, double semiMajorAxis)  // Calculate velocity within an elliptical orbit
        {   // Formula for eliptical orbit velocity. Found at: http://www.braeunig.us/space/orbmech.htm#position (4.45)
            return Math.Sqrt(GravitationalConstant * primaryMass * (2d / DVector3.Distance(primaryPos, bodyPos) - 1d / semiMajorAxis));
        }

        public static double FlightPathAngle(double semiMajorAxis, double eccentricity, Body body)  // angle from cross of position vector and orbit up vector to direction body's travel
        {   // Formula for flight path angle found at: https://www.12000.org/my_notes/dynamics_cheat_sheet/rech3.htm
            double r = (body.Primary.GetWorldPosition() - body.GetWorldPosition()).Magnitude(); // Distance between primary and body
            double fpa = Math.Acos(Math.Sqrt((Math.Pow(semiMajorAxis, 2d) * (1d - Math.Pow(eccentricity, 2d)) / (r * (2d * semiMajorAxis - r)))));
            if (double.IsNaN(fpa)) return 0;
            else return fpa;
        }

        public static Vector3 OrbitVectorUp(double semiMajorAxis, double eccentricity, double inclination, double longitudeOfAscendingNode) // The up vector relative to the orbital plane
        {
            Vector3 trueAnomaly0;   // Direction vector to point in orbit with a true anomaly of 0 degrees
            Vector3 trueAnomaly90;  // 90 degrees

            Quaternion orbitRotation = Quaternion.Euler(Vector3.up * -(float)longitudeOfAscendingNode);
            orbitRotation *= Quaternion.AngleAxis((float)inclination, Vector3.right);   // Calculate rotation of orbit. Argument of periapsis is not needed
                                                                                        // as this parameter rotates the orbit around this axis already
            double d = RadiusAtPointInOrbit(semiMajorAxis, eccentricity, 0);
            double x = Math.Sin(0 * (Math.PI / 180.0d) * (Math.PI / 180.0d));
            double z = Math.Cos(0 * (Math.PI / 180.0d) * (Math.PI / 180.0d));

            trueAnomaly0 = orbitRotation * new Vector3((float)x, 0, (float)z).normalized;   // Record position at true anomaly 0

            d = RadiusAtPointInOrbit(semiMajorAxis, eccentricity, 90);
            x = Math.Sin(90 * (Math.PI / 180.0d) * (Math.PI / 180.0d));
            z = Math.Cos(90 * (Math.PI / 180.0d) * (Math.PI / 180.0d));

            trueAnomaly90 = orbitRotation * new Vector3((float)x, 0, (float)z).normalized;

            return Vector3.Cross(trueAnomaly0, trueAnomaly90).normalized;
        }

        public static Vector3 OrbitVector90(double semiMajorAxis, double eccentricity) // The direction vector of the position in orbit at a true anomaly of 90
        {
            double d = RadiusAtPointInOrbit(semiMajorAxis, eccentricity, 90);
            double x = Math.Sin(90 * Math.PI / 180.0d);
            double z = Math.Cos(90 * Math.PI / 180.0d);

            return new Vector3((float)x, 0, (float)z).normalized;
        }

        public static Vector3 CenterOfEclipse(Orbit orbit, Vector3 primaryGameObjectPos) // Sums the position vectors of apoapsis and periapsis and divides by two
        {
            return (((PositionInOrbit(0d, orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode,
                   orbit.ArgumentOfPeriapsis) / (float)AstronomicalUnit * (float)UniverseManager.GetUniverseScale()) +
                   (PositionInOrbit(180d, orbit.SemiMajorAxis, orbit.Eccentricity, orbit.Inclination,
                   orbit.LongitudeOfAscendingNode, orbit.ArgumentOfPeriapsis) / (float)AstronomicalUnit *
                   (float)UniverseManager.GetUniverseScale())) / 2f) + primaryGameObjectPos;
        }

        public static double HillSphereRadius(DVector3 body1Position, DVector3 body2Position, double body1Mass, double body2Mass)   // Calculates radius of dominant gravitational influence
        {
            double r = DVector3.Distance(body1Position, body2Position);
            double hsr = r * Math.Pow(body1Mass / (3d * body2Mass), 1d / 3d);
            return hsr;
        } 

        public static double PeriapsisRadius(double semiMajorAxis, double eccentricity) // Radius of the periapsis of an orbit
        {
            return semiMajorAxis * (1d - eccentricity);
        }

        public static double ApoapsisRadius(double semiMajorAxis, double eccentricity)  // Radius of the apoapsis of an orbit
        {
            return semiMajorAxis * (1 + eccentricity);
        }

        public static double SemiMinorAxis(double semiMajorAxis, double eccentricity)   // Calculates length of the semi major axis
        {
            return semiMajorAxis * Math.Sqrt(1d - Math.Pow(eccentricity, 2d));
        }
    }   

    public class DVector3   // Double precision 3D vector
    {
        public double x = 0.0d;
        public double y = 0.0d;
        public double z = 0.0d;

        public static DVector3 operator +(DVector3 dv1, DVector3 dv2)
        {
            dv1.x += dv2.x;
            dv1.y += dv2.y;
            dv1.z += dv2.z;

            return dv1;
        }   // Operator Overloads

        public static DVector3 operator -(DVector3 dv1, DVector3 dv2)
        {
            dv1.x -= dv2.x;
            dv1.y -= dv2.y;
            dv1.z -= dv2.z;

            return dv1;
        }

        public static DVector3 operator *(DVector3 dv, double scalar)
        {
            return new DVector3(dv.x * scalar, dv.y * scalar, dv.z * scalar);
        }

        public static DVector3 operator *(double scalar, DVector3 dv)
        {
            return new DVector3(dv.x * scalar, dv.y * scalar, dv.z * scalar);
        }

        public static DVector3 operator /(DVector3 dv, double scalar)
        {
            dv.x /= scalar;
            dv.y /= scalar;
            dv.z /= scalar;

            return dv;
        }

        public double Magnitude()
        {
            return Math.Sqrt((x * x) + (y * y) + (z * z));
        }

        public static double Distance(DVector3 dv1, DVector3 dv2)
        {
            return Math.Sqrt(Math.Pow(dv2.x - dv1.x, 2) + Math.Pow(dv2.y - dv1.y, 2) + Math.Pow(dv2.z - dv1.z, 2));
        }

        public void Normalize() // normalize the vector
        {
            double mag = Magnitude();

            x /= mag;
            y /= mag;
            z /= mag;
        } 

        public DVector3 Normalized()    // Return the vector normalized
        {
            double mag = Magnitude();

            return new DVector3(x / mag, y / mag, z / mag);
        }

        public static double AngleBetweenVectors(DVector3 dv1, DVector3 dv2)
        {
            return Math.Acos((DotProduct(dv1, dv2) / (dv1.Magnitude() * dv2.Magnitude())));
        }

        public static DVector3 CrossProduct(DVector3 dv1, DVector3 dv2)
        {
            return new DVector3((dv1.y * dv2.z) - (dv1.z * dv2.y), (dv1.z * dv2.x) - (dv1.x * dv2.z), (dv1.x * dv2.y) - (dv1.y * dv2.x));
        }

        public static double DotProduct(DVector3 dv1, DVector3 dv2)
        {
            return dv1.x * dv2.x + dv1.y * dv2.y + dv1.z * dv2.z;
        }

        public DVector3(double _x, double _y, double _z)
        {
            x = _x;
            y = _y;
            z = _z;
        }   // Initialise with doubles

        public DVector3(Vector3 v)
        {
            x = (double)v.x;
            y = (double)v.y;
            z = (double)v.z;
        }   // Initialise with float vector3

        public DVector3(DVector3 dv)    // Copy dvector into this vector
        {
            x = dv.x;
            y = dv.y;
            z = dv.z;
        }

        public DVector3()   // Initialise values to 0d
        {
            x = 0d;
            y = 0d;
            z = 0d;
        }

        public static Vector3 GetFloatVector(DVector3 dv)
        {
            return new Vector3((float)dv.x, (float)dv.y, (float)dv.z);
        }

        public static void LogDVector3(DVector3 dv)
        {
            Debug.Log(dv.x + ", " + dv.y + ", " + dv.z);
        }
    } 
}
