module.exports = {
    branches: ['main'],
    plugins: [
        [
            '@semantic-release/commit-analyzer',
            {preset: 'conventionalcommits'}
        ],
        [
            '@semantic-release/release-notes-generator',
            {preset: 'conventionalcommits'},
        ],
        '@semantic-release/changelog',
        [
            '@semantic-release/npm',
            {pkgRoot: 'packages/gsw_cpp'}
        ],
        [
            '@semantic-release/npm',
            {pkgRoot: 'packages/gsw_node'}
        ],
        [
            '@semantic-release/exec',
            {execCwd: 'packages/gsw_cpp', publishCmd: 'pnpm pack'}
        ],
        [
            '@semantic-release/exec',
            {execCwd: 'packages/gsw_node', publishCmd: 'pnpm pack'}
        ],
        [
            '@semantic-release/github',
            {
                assets: [
                    'packages/gsw_cpp/greenroom-robotics-gsw_cpp-*.tgz',
                    'packages/gsw_node/greenroom-robotics-gsw_node-*.tgz'
                ]
            }
        ],
        '@semantic-release/git',
    ],
};
