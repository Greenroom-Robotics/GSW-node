{
	"targets": [
		{
			"includes": [
				"auto.gypi"
			],
			"sources": [
				"TeosBindings.cpp"
			],
			"libraries": [
				"/build/TeosCpp/build/libteos-cpp.a"
			],
			"include_dirs":[
				"../TeosCpp/include"
			]
		}
	],
	"includes": [
		"auto-top.gypi"
	]
}
