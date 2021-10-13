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
				"<(module_root_dir)/TeosCpp/build/libteos-cpp.a"
			],
			"include_dirs":[
				"<(module_root_dir)/TeosCpp/include"
			]
		}
	],
	"includes": [
		"auto-top.gypi"
	]
}
