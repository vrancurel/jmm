
CLASSES = \
MixtureModel/MixtureModel.class

ProbTest.class: ProbTest.java $(CLASSES)
	javac ProbTest.java

$(CLASSES): MixtureModel/*.java
	cd MixtureModel ; javac *.java

clean:
	rm MixtureModel/*.class *.class
