import AssemblyKeys._

name := "phac"

organization := "cmpg"

version := "0.1"

scalaVersion := "2.11.4"

EclipseKeys.withSource := true

exportJars := true

libraryDependencies += "org.slf4j" % "slf4j-api" % "1.6.6"

libraryDependencies += "ch.qos.logback" % "logback-classic" % "1.0.9"

libraryDependencies += "junit" % "junit" % "4.11" % "test"

libraryDependencies += "org.specs2" %%"specs2" % "2.3.13" % "test"
 
libraryDependencies += "org.mockito" % "mockito-core" % "1.9.5" % "test"     
            
libraryDependencies += "net.sf.opencsv" % "opencsv" % "2.3"

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.0"

libraryDependencies += "com.github.scopt" %% "scopt" % "3.3.0"

libraryDependencies += "org.json" % "json" % "20090211"

resolvers += Resolver.sonatypeRepo("public")

mainClass in assembly := Some("be.cmpg.SmallSubnetworkAnalysis")

jarName in assembly := "SSA.jar"

test in assembly := {}
            
