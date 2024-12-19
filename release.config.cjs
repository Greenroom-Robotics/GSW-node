module.exports = {
    branches: ['main'],
    plugins: [
        [
            '@semantic-release/commit-analyzer',
            { preset: 'conventionalcommits' }
        ],
        [
            '@semantic-release/release-notes-generator',
            { preset: 'conventionalcommits' },
        ],
        '@semantic-release/changelog',
        [
            '@semantic-release/npm',
            { pkgRoot: 'packages/gsw_cpp' }
        ],
        [
            '@semantic-release/npm',
            { pkgRoot: 'packages/gsw_node' }
        ],
        '@semantic-release/github',
        '@semantic-release/git',
    ],
};
