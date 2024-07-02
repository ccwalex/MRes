run("Z Project...", "projection=[Sum Slices]");
main = getTitle();
run("Arrange Channels...");
selectImage(main);

//run("Channels Tool...");

Property.set("CompositeProjection", "Sum");
Stack.setDisplayMode("composite");
Stack.setActiveChannels("011");

run("Scale Bar...", "width=100 height=20 thickness=10 font=36 bold overlay");
run("Stack to RGB");


selectImage(main)
Stack.setActiveChannels("110");

run("Scale Bar...", "width=100 height=20 thickness=10 font=36 bold overlay");
run("Stack to RGB");

selectImage(main)
Stack.setActiveChannels("101");

run("Scale Bar...", "width=100 height=20 thickness=10 font=36 bold overlay");
run("Stack to RGB");
selectImage(main);
Stack.setActiveChannels("111");

run("Scale Bar...", "width=100 height=20 thickness=10 font=36 bold overlay");
run("Stack to RGB");
