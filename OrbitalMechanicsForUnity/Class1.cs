using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UnityEngine;
using UnityEditor;

namespace OrbitalMechanicsForUnity
{
    [Serializable]
    public class Body : MonoBehaviour
    { 
        public double mass = 5972000000000000000000000d;
        public bool DrawHillSphere = true;
        public double hillSphereRadius;
        public Body Primary = null;

        [SerializeField]
        Orbit orbit = new Orbit();
        DVector3 position;

        public void UpdateBody()
        {
            if(UniverseManager.GetTimeScale() != 0 && Primary)
            {

                //double av = OrbitalMechanics.AngularVelocity(Primary.mass, Primary.GetPosition(),
                //                                            position) *
                //                                            (180d / Math.PI) * UniverseManager.GetTimeScale() * Time.deltaTime;

                orbit.TrueAnomaly = OrbitalMechanics.TrueAnomalyAfterTime(orbit.SemiMajorAxis, orbit.Eccentricity, orbit.TrueAnomaly,
                                                                          Time.deltaTime * (float)UniverseManager.GetTimeScale(), Primary.mass);

                //orbit.TrueAnomaly += av;

                if(orbit.TrueAnomaly > 360) // Clamp the true anomaly 0-360 to prevent the value from getting too large and causing data issues
                {
                    orbit.TrueAnomaly -= 360;
                }
                else if (orbit.TrueAnomaly < 0)
                {
                    orbit.TrueAnomaly += 360;
                }

                position = GetWorldPosition();
            }
        }
        private void Start()
        {
            position = GetWorldPosition();

            if(UniverseManager.scaled && Primary)
            {
                gameObject.transform.localPosition = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit)) / (float)OrbitalMechanics.AstronomicalUnit
                                                   * (float)UniverseManager.GetUniverseScale();
            }
        }
        private void Update()
        {
            UpdateBody();

            if (UniverseManager.scaled && Primary)
            {
                gameObject.transform.localPosition = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit)) / (float)OrbitalMechanics.AstronomicalUnit
                                                   * (float)UniverseManager.GetUniverseScale();
            }
        }
        public DVector3 GetPosition()
        {
            return position;
        }

        public DVector3 GetWorldPosition()
        {
            if(Primary)
            {
                //return OrbitalMechanics.PositionInOrbit(orbit) + (new DVector3(Primary.transform.position) * OrbitalMechanics.AstronomicalUnit);
                if(Primary.GetPosition() == null)
                {
                    //Debug.Log("Primary Position is null for " + name);
                    return OrbitalMechanics.PositionInOrbit(orbit) + OrbitalMechanics.PositionInOrbit(Primary.Orbit());
                }
                else
                {
                    return OrbitalMechanics.PositionInOrbit(orbit) + Primary.GetPosition();
                }
                
            }
            else
            {
                return new DVector3(transform.position) * OrbitalMechanics.AstronomicalUnit;
            }
        }

        public Orbit Orbit()
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
    }

    [CustomEditor(typeof(Body))]    // informs Unity which component it should act as an editor for
    [CanEditMultipleObjects]        // Tells Unity that multiple objects can be selected and edited at the same time
    public class BodyEditor : Editor    // Editor for body class
    {
        SerializedProperty mass, orbit, drawHillSphere, primary;

        void OnEnable()
        {
            mass = serializedObject.FindProperty("mass");
            orbit = serializedObject.FindProperty("orbit");
            drawHillSphere = serializedObject.FindProperty("DrawHillSphere");
            primary = serializedObject.FindProperty("Primary");
        }

        public override void OnInspectorGUI()
        {
            serializedObject.Update();

            EditorGUILayout.ObjectField(primary);
            EditorGUILayout.PropertyField(mass);

            Body myBody = (Body)target;

            if(myBody.Primary)
            {
                EditorGUILayout.PropertyField(drawHillSphere);
                EditorGUILayout.PropertyField(orbit);
            }

            serializedObject.ApplyModifiedProperties();
        }

        public void OnSceneGUI()
        { 
            var body = target as Body;

            if (!body) return;          

            if(body.Primary)
            {
                Orbit orbit = body.Orbit();
                
                body.hillSphereRadius = OrbitalMechanics.HillSphereRadius(new DVector3(body.transform.position) * OrbitalMechanics.AstronomicalUnit, // Calculate hill sphere radius
                                                                          new DVector3(body.Primary.transform.position) * OrbitalMechanics.AstronomicalUnit, body.mass, body.Primary.mass);

                if (orbit.RenderOrbit)   // Render the orbit path
                {
                    var transform = body.gameObject.transform;

                    Vector3 offset = new Vector3(0, 0, 0);  // Vector that states where to draw the orbit relative to.
                    
                    offset = body.Primary.transform.position;
                    
                    Vector3 v1 = OrbitalMechanics.PositionInOrbit(0d, orbit.SemiMajorAxis,
                                 orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode,
                                 orbit.ArgumentOfPeriapsis) / (float)OrbitalMechanics.AstronomicalUnit * (float)UniverseManager.GetUniverseScale() + offset;

                    if (!Application.isPlaying)
                    {
                        transform.localPosition = DVector3.GetFloatVector(OrbitalMechanics.PositionInOrbit(orbit)) / (float)OrbitalMechanics.AstronomicalUnit
                                                                 * (float)UniverseManager.GetUniverseScale();
                    }

                    Vector3 v2;

                    for (int i = 1; i < 30; i++)
                    {
                        v2 = OrbitalMechanics.PositionInOrbit(((360d / 30d) * i), orbit.SemiMajorAxis,
                                                     orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode,
                                                     orbit.ArgumentOfPeriapsis) / (float)OrbitalMechanics.AstronomicalUnit
                                                     * (float)UniverseManager.GetUniverseScale() + offset;

                        Handles.DrawLine(v1, v2);

                        v1 = v2;
                    }

                    Handles.DrawLine(v1, OrbitalMechanics.PositionInOrbit(0d, orbit.SemiMajorAxis,
                                         orbit.Eccentricity, orbit.Inclination, orbit.LongitudeOfAscendingNode,
                                         orbit.ArgumentOfPeriapsis) / (float)OrbitalMechanics.AstronomicalUnit
                                         * (float)UniverseManager.GetUniverseScale() + offset);
                }
            }
        }
    }

    [Serializable]
    public class Orbit
    {
        public double SemiMajorAxis;
        [Range(0f, 1f)]
        public double Eccentricity;
        [Range(-180f, 180f)]
        public double Inclination;
        [Range(0f, 360f)]
        public double LongitudeOfAscendingNode;
        [Range(0f, 360f)]
        public double ArgumentOfPeriapsis;
        [Range(0f, 360f)]
        public double TrueAnomaly;

        public bool RenderOrbit;

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
    }

    // Orbit editor class??

    static public class UniverseManager
    {
        static double TimeScale = 1d;
        static double UniverseScale = 10d;
        public static bool scaled = false;
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
    }

    public class OrbitalMechanics
    {
        public const double GravitationalConstant = 0.0000000000667408d;
        static public double AstronomicalUnit = 149597870700d / UniverseManager.GetUniverseScale();
        static public double SunMass = 1989000000000000000000000000000d; // TEMP VARIABLE FOR TESTING

        public static double RadiusAtPointInOrbit(double semiMajorAxis, double eccentricity, double trueAnomaly)
        {
            // Debug.Log("(" + semiMajorAxis + " * (1 - (e^2))) / (1 + (" + eccentricity + " * cos(" + trueAnomalyInRad + ")");
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

        public static double AngularVelocity(double primaryMass, DVector3 primaryPosition, DVector3 bodyPosition) // Calculate a bodies angular belocity in radians.
        {
            double r = DVector3.Distance(primaryPosition, bodyPosition);
            return Math.Sqrt((GravitationalConstant * primaryMass) / Math.Pow(r, 3d));
        }

        public static double TrueAnomalyAfterTime(double semiMajorAxis, double eccentricity, double trueAnomaly, float timePassed, double primaryMass)
        {
            double v0 = trueAnomaly * Math.PI / 180d;   // Initial true anomaly in radians
            Debug.Log("v0: " + v0);
            double e0 = Math.Acos((eccentricity + Math.Cos(v0)) / (1 + eccentricity * Math.Cos(v0))); // Initial eccentric anomaly
            Debug.Log("e0: " + e0);
            double m0 = e0 - eccentricity * Math.Sin(e0);   // Initial mean anomaly
            Debug.Log("m0: " + m0);
            double n = Math.Sqrt((GravitationalConstant * primaryMass) / Math.Pow(semiMajorAxis, 3));   // Mean motion
            Debug.Log("n: " + n);
            double m = m0 + n * (timePassed);   // mean anomaly at new position
            Debug.Log("m: " + m);
            double a1 = m + eccentricity * Math.Sin(2d);    // first answer to iteration
            double a2 = 0; 
            bool iterationComplete = false;
            while(!iterationComplete)
            {
                a2 = m + eccentricity * Math.Sin(a1);

                if(Math.Round(a1, 5) == Math.Round(a2, 5))
                {
                    iterationComplete = true;
                    Debug.Log("E: " + a2);
                }
                else
                {
                    a1 = a2;
                }
            }
            
            double v = Math.Acos((Math.Cos(a2) - eccentricity) / (1 - eccentricity * Math.Cos(a2)));
            Debug.Log("v: " + v);
            return v * 180d / Math.PI;

        }

        public static Vector3 PositionInOrbit(double trueAnomaly, double semiMajorAxis, double eccentricity, double inclination, double longitudeOfAscendingNode,
                                              double argumentOfPeriapsis) // returns vector3 of position in orbit at given true anomaly
        {   // CREATE DOUBLE PRECISION QUATERNION CLASS?
            double a = trueAnomaly * Math.PI / 180.0d;
            double d = RadiusAtPointInOrbit(semiMajorAxis, eccentricity, a);
            double x = Math.Sin(trueAnomaly * Math.PI / 180.0d);
            double z = Math.Cos(trueAnomaly * Math.PI / 180.0d);

            Quaternion pointRotation = Quaternion.Euler((float)longitudeOfAscendingNode * Vector3.up);
            pointRotation *= Quaternion.AngleAxis((float)inclination, OrbitVector90(semiMajorAxis, eccentricity));
            pointRotation *= Quaternion.AngleAxis((float)-argumentOfPeriapsis, Vector3.up);

            return pointRotation * new Vector3((float)x, 0, (float)z).normalized * (float)(d);
        }

        public static DVector3 PositionInOrbit(Orbit orbit) // returns DVector3 of position in orbit
        {   // CREATE DOUBLE PRECISION QUATERNION CLASS?
            double a = orbit.TrueAnomaly * Math.PI / 180.0d;
            double d = RadiusAtPointInOrbit(orbit.SemiMajorAxis, orbit.Eccentricity, a);
            double x = Math.Sin(orbit.TrueAnomaly * Math.PI / 180.0d);
            double z = Math.Cos(orbit.TrueAnomaly * Math.PI / 180.0d);

            Quaternion pointRotation = Quaternion.Euler((float)orbit.LongitudeOfAscendingNode * Vector3.up);
            pointRotation *= Quaternion.AngleAxis((float)orbit.Inclination, OrbitVector90(orbit.SemiMajorAxis, orbit.Eccentricity));
            pointRotation *= Quaternion.AngleAxis((float)-orbit.ArgumentOfPeriapsis, Vector3.up);

            return new DVector3(pointRotation * new Vector3((float)x, 0, (float)z).normalized * (float)(d));
        }

        public static Vector3 OrbitVectorUp(double semiMajorAxis, double eccentricity, double inclination, double longitudeOfAscendingNode) // The up vector relative to the orbital plane
        {
            Vector3 trueAnomaly0;   // Direction vector to point in orbit with a true anomaly of 0 degrees
            Vector3 trueAnomaly90;  // 90 degrees

            Quaternion orbitRotation = Quaternion.Euler(Vector3.up * (float)longitudeOfAscendingNode);
            orbitRotation *= Quaternion.AngleAxis((float)inclination, Vector3.right);

            double d = RadiusAtPointInOrbit(semiMajorAxis, eccentricity, 0);
            double x = Math.Sin(0 * (Math.PI / 180.0d) * (Math.PI / 180.0d));
            double z = Math.Cos(0 * (Math.PI / 180.0d) * (Math.PI / 180.0d));

            trueAnomaly0 = orbitRotation * new Vector3((float)x, 0, (float)z).normalized;

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

        public static double HillSphereRadius(DVector3 body1Position, DVector3 body2Position, double body1Mass, double body2Mass)
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
    }

    public class DVector3
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

        public void Normalize()
        {
            double mag = Magnitude();

            x /= mag;
            y /= mag;
            z /= mag;
        }

        public DVector3 Normalized()
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
